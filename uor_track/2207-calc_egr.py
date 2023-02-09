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

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

path = '/home/ys17-23/Extension2/renql/ERA5_mon'
figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig/"
outdir = "/home/ys17-23/Extension2/renql/project/uor_track/mdata"
#lonl=0 
#lonr=359
#lats=0
#latn=90
lonl=20 
lonr=140
lats=15
latn=70
ilat = np.arange(lats, latn+0.1, 1)
ilon = np.arange(lonl, lonr+0.1, 1)
lat_sp = 20
lon_sp = 30 #60 #
ga = 9.80665 # Gravitational acceleration
a  = 6378388 # the radius of earth, m
numod= [chr(i) for i in range(97,115)]

def main_run():
    #varname = 'th';dvar=varname;scale=1;unit='K';cnlev=np.arange(215,295.1,5)
    varname = 'dthdy';dvar=varname;scale=-100000;unit='K/100km';cnlev=np.arange(-1.6,1.61,0.2)
    #varname = 'dvdz';dvar=varname;scale=1000;unit='$10^3 s^{-1}$';cnlev=np.arange(0,8.1,0.5)
    #varname = 'egr';dvar='EGR';scale=1;unit='$day^{-1}$';cnlev=np.arange(0,1.61,0.1)
    #varname = 'eff_egr';dvar='EGR_eff';scale=1;unit='$day^{-1}$';cnlev=np.arange(0,1.61,0.1)
    #varname = 'N2';dvar=varname;scale=100000;unit='10$^{-5}$ s$^{-2}$';cnlev=np.arange(1,34,2)
    #varname = 'dtdz';dvar=varname;scale=1000;unit='K/km';cnlev=np.arange(1,9.1,0.5)
    #varname = 'eff_dtdz';dvar=varname;scale=1000;unit='K/km';cnlev=np.arange(1,9.1,0.5)
    #varname = 'dtdz_moist';dvar=varname;scale=1000;unit='K/km';cnlev=np.arange(8.5,11.8,0.2)
    outfile = '%s/month41_%s.nc'%(outdir,varname)
    #calc_monthly_egr(outfile,varname)
    #calc_monthly_n2(outfile,varname)
    calc_monthly_dtdy(outfile,varname)
    #alc_monthly_dudz(outfile,varname)
    draw_season_4x3(outfile,varname,scale,unit,cnlev,dvar)

def calc_monthly_dtdy(outfile,varname):
    level = [250,500,850]
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset("%s/u/ERA5_u_1980.nc"%path)
    var = ds['u'].sel(level=level,latitude=ilat,longitude=ilon).load()
    var = var.groupby(var.time.dt.month).mean('time')
    var.data = np.zeros(var.data.shape)
    del ds
    f0  = 2*(2*np.pi/24.0/3600.0)*np.sin(ilat*np.pi/180.0)
    
    nyear = 0
    for year in range(1980,2021,1):
        print('calc for %d'%year)
        ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon
                ).groupby(ds.time.dt.month).mean('time').data
        if varname == ['th','dthdy']:
            t = dynamic_calc.calc_pot_temp(t, np.array(level), 1)
        if varname in ['teqv','dteqvdy']:
            ds = xr.open_dataset("%s/q/ERA5_q_%d.nc"%(path,year))
            q = ds['q'].sel(level=level,latitude=ilat,longitude=ilon
                    ).groupby(ds.time.dt.month).mean('time').data
            t = dynamic_calc.calc_teqv(t,q)
            del q
            t = dynamic_calc.calc_pot_temp(t, np.array(level), 1)
        del ds
        gc.collect()
        if varname in ['dtdy','dthdy','dteqvdy']:
            var.data = var.data + dynamic_calc.center_diff(t, ilat*np.pi/180.0, 2)/a
        if varname in ['t','th','teqv']:
            var.data = var.data + t
        nyear = nyear + 1
 
    print('total year %d'%nyear)
    var.data = var.data/nyear
    var.attrs['units'] = 'K/m'
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def calc_monthly_dudz(outfile,varname):
    #level = [225, 250, 300, 450, 500, 550, 825, 850, 875]
    level = [225, 250, 450, 500, 825, 850]
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset("%s/u/ERA5_u_1980.nc"%path)
    var = ds['u'].sel(level=level,latitude=ilat,longitude=ilon).load()
    var = var.groupby(var.time.dt.month).mean('time')
    var.data = np.zeros(var.data.shape)
    del ds
    
    nyear = 0
    for year in range(1980,2021,1):
        print('calc for %d'%year)
        ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon
            ).groupby(ds.time.dt.dayofyear).mean('time').data/ga
        
        ds = xr.open_dataset("%s/v/ERA5_v_%d.nc"%(path,year))
        u = ds['v'].sel(level=level,latitude=ilat,longitude=ilon
                ).groupby(ds.time.dt.dayofyear).mean('time'
                ).rename({'dayofyear':'time'})
        u.coords['time'] = pd.date_range(start='%d-01-01'%year,
                end='%d-12-31'%year,freq='1D',closed=None)
        
        if varname == 'duvdz':
            ds = xr.open_dataset("%s/v/ERA5_v_%d.nc"%(path,year))
            v = ds['v'].sel(level=level,latitude=ilat,longitude=ilon
                ).groupby(ds.time.dt.dayofyear).mean('time').data
            del ds
            #u.data = np.sqrt(np.power(u.data,2)+np.power(v,2))
            u.data = np.sqrt(np.power(dynamic_calc.upward_diff_z(u.data, z),2)+
                    np.power(dynamic_calc.upward_diff_z(v, z),2))
            del v
        else:
            del ds
            u.data = np.abs(dynamic_calc.upward_diff_z(u.data, z))
        gc.collect()
        u = u.groupby(u.time.dt.month).mean('time')
        var.data = var.data + u.data 
        nyear = nyear + 1
 
    print('total year %d'%nyear)
    var.data = var.data/nyear
    var.attrs['units'] = '1/day'
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def calc_monthly_n2(outfile,varname):
    level = [225, 250, 450, 500, 825, 850]
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset("%s/u/ERA5_u_1980.nc"%path)
    var = ds['u'].sel(level=level,latitude=ilat,longitude=ilon).load()
    var = var.groupby(var.time.dt.month).mean('time')
    var.data = np.zeros(var.data.shape)
    del ds
    nyear = 0
    for year in range(1980,2021,1):
        print('calc for %d'%year)
        ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon
            ).data/ga
        #z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon
        #    ).groupby(ds.time.dt.dayofyear).mean('time').data/ga
        
        ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon)
        #t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon
        #        ).groupby(ds.time.dt.dayofyear).mean('time'
        #        ).rename({'dayofyear':'time'})
        #t.coords['time'] = pd.date_range(start='%d-01-01'%year,
        #        end='%d-12-31'%year,freq='1D',closed=None)
        del ds
        gc.collect()
        
        if varname in ['eff_dtdz',]:
            ds = xr.open_dataset("%s/w/ERA5_w_%d.nc"%(path,year))
            w = ds['w'].sel(level=level,latitude=ilat,longitude=ilon
                ).data
            #w = ds['w'].sel(level=level,latitude=ilat,longitude=ilon
            #    ).groupby(ds.time.dt.dayofyear).mean('time').data
            t.data = dynamic_calc.calc_eff_dzdt(t.data,z,level,w,False)
        else:
            t.data = dynamic_calc.calc_n2(t.data,z,level,varname,False)
        
        t = t.groupby(t.time.dt.month).mean('time')
        var.data = var.data + t.data 
        nyear = nyear + 1
 
    print('total year %d'%nyear)
    var.data = var.data/nyear
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def calc_monthly_egr(outfile,varname):
    #level = [225, 250, 300, 450, 500, 550, 825, 850, 875]
    level = [225, 250, 450, 500, 825, 850]
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    ds = xr.open_dataset("%s/u/ERA5_u_1980.nc"%path)
    var = ds['u'].sel(level=level,latitude=ilat,longitude=ilon).load()
    var = var.groupby(var.time.dt.month).mean('time')
    var.data = np.zeros(var.data.shape)
    del ds
    f0  = 2*(2*np.pi/24.0/3600.0)*np.sin(ilat*np.pi/180.0)
    nyear = 0
    for year in range(1980,2021,1):
        print('calc for %d'%year)
        ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon
            ).groupby(ds.time.dt.month).mean('time').data/ga
        ds = xr.open_dataset("%s/u/ERA5_u_%d.nc"%(path,year))
        u = ds['u'].sel(level=level,latitude=ilat,longitude=ilon
                ).groupby(ds.time.dt.month).mean('time')
        ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon
            ).groupby(ds.time.dt.month).mean('time').data
        
        #ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        #z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon
        #    ).groupby(ds.time.dt.dayofyear).mean('time').data/ga
        #
        #ds = xr.open_dataset("%s/u/ERA5_u_%d.nc"%(path,year))
        #u = ds['u'].sel(level=level,latitude=ilat,longitude=ilon
        #        ).groupby(ds.time.dt.dayofyear).mean('time'
        #        ).rename({'dayofyear':'time'})
        #u.coords['time'] = pd.date_range(start='%d-01-01'%year,
        #        end='%d-12-31'%year,freq='1D',closed=None)
        #
        #ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        #t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon
        #    ).groupby(ds.time.dt.dayofyear).mean('time').data
        
        #ds = xr.open_dataset("%s/v/ERA5_v_%d.nc"%(path,year))
        #v = ds['v'].sel(level=level,latitude=ilat,longitude=ilon
        #    ).groupby(ds.time.dt.dayofyear).mean('time').data
        
        #ds = xr.open_dataset("%s/t/ERA5_t_%d.nc"%(path,year))
        #t = ds['t'].sel(level=level,latitude=ilat,longitude=ilon).data
        #ds = xr.open_dataset("%s/z/ERA5_z_%d.nc"%(path,year))
        #z = ds['z'].sel(level=level,latitude=ilat,longitude=ilon).data/ga
        #ds = xr.open_dataset("%s/u/ERA5_u_%d.nc"%(path,year))
        #u = ds['u'].sel(level=level,latitude=ilat,longitude=ilon)
        #ds = xr.open_dataset("%s/v/ERA5_v_%d.nc"%(path,year))
        #v = ds['v'].sel(level=level,latitude=ilat,longitude=ilon)
        del ds
        gc.collect()
        if varname in ['eff_egr',]:
            ds = xr.open_dataset("%s/w/ERA5_w_%d.nc"%(path,year))
            w = ds['w'].sel(level=level,latitude=ilat,longitude=ilon).data
            u.data = dynamic_calc.calc_eff_egr(t, z, u.data, f0, level, w)
            #u.data = dynamic_calc.calc_eff_egr(t, z, u.data, f0, level, w, v)
        else:
            u.data = dynamic_calc.calc_egr(t, z, u.data, f0, level)
        #u = u.groupby(u.time.dt.month).mean('time')
        var.data = var.data + u.data 
        nyear = nyear + 1
 
    print('total year %d'%nyear)
    var.data = var.data*86400/nyear
    var.attrs['long_name'] = 'Eady growth rate'
    var.attrs['units'] = '1/day'
    print(var)
    ds1 = var.to_dataset(name=varname)
    ds1.to_netcdf(outfile,'w')

def draw_season_4x3(outfile,varname,scal,unit,cnlev,dvar):
    lev = [850,500,250]
    titls= ['DJF','MAM','JJA','SON']
    
    ds = xr.open_dataset(outfile)
    #lat = ds.latitude
    #lon = ds.longitude
    #ilon = lon[(lon>=lonl) & (lon<=lonr)]
    #ilat = lat[(lat>=lats) & (lat<=latn)]
    var = scal*ds[varname].sel(level=lev,longitude=ilon,latitude=ilat).data
    '''
    ds = xr.open_dataset('%s/ERA-Interim_EGR-year.nc'%outdir)
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var = 24*3600*ds['month_ave'].sel(lev=lev,lon=ilon,lat=ilat).data.mean(axis=0)
    ''' 

    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.3 #0.25 #
    
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
    
    ds = xr.open_dataset('%s/ERA5_mon_u_1979-2020.nc'%path)
    da = ds['u'].sel(level=200,longitude=ilon,
        latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").data
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    gc.collect()

    for nl in range(0,3,1):
        for nm in range(0,nrow,1):
            if nm == 0:
                shad = (var[0,nl,:,:]+var[1,nl,:,:]+var[11,nl,:,:])/3.0
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                shad = np.mean(var[(3*nm-1):(3*nm+2),nl,:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
            print('%d %s %s : min = %f ; max = %f'%(lev[nl],varname,titls[nm],
                np.nanmin(shad),np.nanmax(shad)))
            if nl==0:
                shad = np.ma.array(shad, mask=(phis>1500))
            axe = ax[nm][nl]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
                edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title("(%s) %dhPa %s"%(numod[3*nm+nl],lev[nl],titls[nm]),
                fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, shad, cnlev, 
                 transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
                 transform=ccrs.PlateCarree(),colors=jetcolor,linewidths=2.2)
            
            if nl == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.35, bmlo+0.005, 0.55, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.15,bmlo-0.005, '%s (%s)'%(dvar,unit),fontsize=title_font,
        horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig('%s/seasonal_%s.png'%(figdir,varname), bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

