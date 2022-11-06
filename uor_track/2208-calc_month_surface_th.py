#!/usr/bin/env python
import sys
import subprocess, os
import xarray as xr
import numpy as np
import pandas as pd 
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from renql import dynamic_calc 
#matplotlib.use('Agg')

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon'
outdir = '/home/users/qd201969/uor_track/mdata'
figdir = '/home/users/qd201969/uor_track/fig'
lonl=15 
lonr=145
lats=15
latn=70
ilat = np.arange(lats, latn+0.1, 0.25)
ilon = np.arange(lonl, lonr+0.1, 0.25)
lat_sp = 20
lon_sp = 30 #60 #
ga = 9.80665 # Gravitational acceleration
a  = 6378388 # the radius of earth, m

def main_run():
    #varname = '2m_th';dvar=varname;scale=1;unit='K';cnlev=np.arange(245,325.1,5)
    varname = 'd2mthdy';dvar=varname;scale=-100000;unit='K/100km';cnlev=np.arange(-1.6,1.61,0.2)
    outfile = '%s/month41_%s.nc'%(outdir,varname)
    #calc_monthly_surface(outfile,varname)
    calc_monthly_dtdy(outfile,varname)
    draw_season_4x1(outfile,varname,scale,unit,cnlev,dvar)

def calc_monthly_dtdy(outfile,varname):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset('%s/month41_2m_th.nc'%outdir)
    var = ds['2m_th']
    del ds
    gc.collect()
    
    var.data = dynamic_calc.center_diff(var.data, 
            ilat*np.pi/180.0, 1)/a
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def calc_monthly_surface(outfile,varname):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset('%s/ERA5_mon_t2m_1979-2020.nc'%path)
    var = ds['t2m'].sel(longitude=ilon,latitude=ilat).groupby(
            ds.time.dt.month).mean('time')
    ds = xr.open_dataset('%s/ERA5_mon_sp_1979-2020.nc'%path)
    sp = ds['sp'].sel(longitude=ilon,latitude=ilat).groupby(
            ds.time.dt.month).mean('time').data
    del ds
    gc.collect()
    var.data = var.data*np.power(100000.0/sp,0.286)
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def draw_season_4x1(outfile,varname,scal,unit,cnlev,dvar):
    titls= ['DJF','MAM','JJA','SON']
    
    ds = xr.open_dataset(outfile)
    var = scal*ds[varname].sel(longitude=ilon,latitude=ilat).data

    nrow = 4 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
   
    if cnlev[0] >=0 :
        #ncmap = colors.ListedColormap(cmaps.topo_15lev(range(0,16,1))[::-1])
        ncmap = cmaps.precip2_17lev
    else:
        ncmap = cmaps.BlueDarkRed18  
    norm = colors.BoundaryNorm(boundaries=cnlev,
        ncolors=ncmap.N,extend='both')
    jetcolor = 'darkviolet'
    
    uwndpath = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
    ds = xr.open_dataset(uwndpath)
    da = ds['u'].sel(level=200,longitude=ilon,
        latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").data
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    gc.collect()

    for nm in range(0,nrow,1):
        if nm == 0:
            shad = (var[0,:,:]+var[1,:,:]+var[11,:,:])/3.0
            uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
        else:
            shad = np.mean(var[(3*nm-1):(3*nm+2),:,:],axis=0)
            uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
        print('%s %s : min = %f ; max = %f'%(varname,titls[nm],
            np.nanmin(shad),np.nanmax(shad)))
        axe = ax[nm]
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
            edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title("%s %s"%(titls[nm],dvar
            ),fontsize=title_font,fontdict=font)

        cont = axe.contourf(ilon, ilat, shad, cnlev, 
             transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000], 
             transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        #jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
        #     transform=ccrs.PlateCarree(),colors=jetcolor,linewidths=2.2)
        
        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nm == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.65, bmlo+0.05, 0.01, 1-bmlo-0.1]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
    plt.tight_layout(w_pad=0.1,rect=(0,bmlo,1,1))
    plt.savefig('%s/seasonal_%s.png'%(figdir,varname), bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

