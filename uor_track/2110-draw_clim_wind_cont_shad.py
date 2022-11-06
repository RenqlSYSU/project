#!/usr/bin/env python
'''
read uwnd, vwnd, z to draw monthly wind (vector) and geopotential (shaded)
Loop through the height (850, 500, 250)
maybe later the shaded variable can be changed for t, PV, dtdy

20211007
'''
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from renql import dynamic_calc

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30
lev = [850,500,250]

nrow = 4
ncol = 3
bmlo = 0.4
BIGFONT=22
MIDFONT=14
SMFONT=10

titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
varname = ['z']
drawvar = ['z']
unit    = ['m']
vcref =[10,20,30] # different levels 
cnlvl =3*[[-4,0.5]]
#cnlvl =[[1300 ,20 ],
#        [5050 ,60 ],
#        [9300 ,100]]
q_mis=15
figdir = "/home/users/qd201969/uor_track/fig/"
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon'

f = xr.open_dataset('%s/ERA5_mon_z_1979-2020.nc'%path)
lat = f.latitude.data
lon = f.longitude.data
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m

for nl in range(0,len(lev),1):
    #da = f['z'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    #var = da.groupby(da.time.dt.month).mean('time')
    #var.data = var.data/9.8
    #print(var)

    ds = xr.open_dataset('%s/ERA5_mon_u_1979-2020.nc'%path)
    da = ds['u'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data

    ds = xr.open_dataset('%s/ERA5_mon_v_1979-2020.nc'%path)
    da = ds['v'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
    vwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da
    gc.collect()
    
    var = dynamic_calc.calc_uv2vr_cfd(uwnd,vwnd,ilat,ilon) 
    #speed = np.power((np.power(uwnd.data,2)+np.power(vwnd.data,2)),0.5)
    #vcref = int(np.percentile(speed,50))
    #maxlvl = int(np.percentile(var.data,80))
    #minlvl = int(np.percentile(var.data,20))
    #if minlvl < 0 :
    if cnlvl[nl][0] < 0 :
        fcolors = cmaps.BlueDarkRed18
    else:
        fcolors = cmaps.precip2_17lev
    #spacig = round((maxlvl-minlvl)/fcolors.N)
    #cnlevels = np.arange(minlvl,maxlvl,spacig)
    cnlevels = np.arange(cnlvl[nl][0], cnlvl[nl][0]+cnlvl[nl][1]*(fcolors.N-1), cnlvl[nl][1])
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

    fig = plt.figure(figsize=(12,12),dpi=100)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
    nm = -1
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            nm = nm+1 
            axe = ax[nr][nc]
            #axe.add_feature(cfeat.LAND.with_scale('110m'), edgecolor='black', linewidth=0.8, zorder=1) 
            axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            #axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=2)
            axe.set_title(str(lev[nl])+" "+titls[nm],fontsize=SMFONT)

            print('min:%f ; max:%f'%(np.nanmin(var[nm,:,:]),np.nanmax(var[nm,:,:])))
            shad = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels,
                         transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
            
            wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd[nm,::q_mis,::q_mis],vwnd[nm,::q_mis,::q_mis],
                    pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="dimgray",
                    width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

            #cont = axe.contour(ilon, ilat, var[nm,:,:], np.arange(1000,15000,cnlvl[nl][1]), 
            #             transform=ccrs.PlateCarree(), colors='darkviolet', linewidths=2)

            topo = axe.contour(ilon, ilat, phis, [1500,3000],
                         transform=ccrs.PlateCarree(),colors='black',linewidths=1.2)

            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.18, bmlo, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)
    axe.quiverkey(wind, 0.92, bmlo-0.01, vcref[nl], r'$%d m/s$'%vcref[nl], labelpos='N',coordinates='figure')

    plt.figtext(0.02,bmlo-0.005, "%dhPa %s (%s)"%(lev[nl],drawvar[0],unit[0]), fontsize=MIDFONT,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir+"%s_%d.png"%(drawvar[0],lev[nl]), bbox_inches='tight',pad_inches=0.01)

