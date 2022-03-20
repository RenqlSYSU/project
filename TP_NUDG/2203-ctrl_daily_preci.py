import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

filein = '/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.'
fileout= '/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/CTRL-Clim_daily_preci.nc'
figdir = '/home/ys17-19/renql/project/TP_NUDG/analysis/fig/CTRL-Clim_daily_preci'

lats = 15 
latn = 55 
lonl = 50  
lonr = 130

def main_run():
    draw_map()

def draw_map():
    lat_sp = 20
    lon_sp = 30 #60 #
    title_font=14
    label_font=10
    
    ds = xr.open_dataset(fileout)
    ilon = ds.lon
    ilat = ds.lat
    var = ds['preci'].resample(time="5D").mean()
    var = var*3600*24*1000 #mm/day
    time = var.indexes['time'].to_datetimeindex()

    ds = xr.open_dataset("/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    
    cnlevels = np.arange(1,18,1)
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    for nt in range(len(time)):
        title = 'CESM AMIP %02d-%02d preci (mm/day)'%(time[nt].month,time[nt].day)
        
        fig = plt.figure(figsize=(12,9),dpi=100)
        axe = plt.axes(projection=ccrs.PlateCarree(central_longitude=180.0))
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
            edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title(title,fontsize=title_font)

        cont = axe.contourf(var.lon, var.lat, var[nt,:,:], cnlevels, 
            transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        topo = axe.contour(phis.lon, phis.lat, phis, [1500,3000], 
            transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)

        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
        axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
        
        cb = plt.colorbar(cont)
        plt.savefig('%s%d'%(figdir,nt), bbox_inches='tight',pad_inches=0.01)

def calc_clim_daily():
    ds = xr.open_dataset('%s1979.daily.PRECL.nc'%filein)
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var = ds['PRECL'].sel(lon=ilon,lat=ilat)
    ds = xr.open_dataset('%s1979.daily.PRECC.nc'%filein)
    var = var + ds['PRECC'].sel(lon=ilon,lat=ilat).data

    for ny in range(1980,2006):
        print(ny)
        ds = xr.open_dataset('%s%d.daily.PRECC.nc'%(filein,ny))
        var = var + ds['PRECC'].sel(lon=ilon,lat=ilat).data
        ds = xr.open_dataset('%s%d.daily.PRECL.nc'%(filein,ny))
        var = var + ds['PRECL'].sel(lon=ilon,lat=ilat).data
    var = var/27
    print(var)
    ds2 = var.to_dataset(name='preci')
    ds2.to_netcdf(fileout,'w')

if __name__=='__main__':
    main_run()
