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

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

lonl=45  #0  #
lonr=135#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30
lev = [850,500,250]
titls= ['DJF','MAM','JJA','SON']
numod= [chr(i) for i in range(97,115)]
vcref =[10,20,30] # different levels 
q_mis=15

figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig/"
path = '/home/ys17-23/Extension2/renql/ERA5_mon'
def main_run():
    cnlvl2=[20,60,100] # contour
    #read_draw_seasonal_4x3('vo',100000,[-3.5,0.5],'vor','s-1')
    #read_draw_seasonal_4x3('w',1,[-0.08,0.01],'omega','Pa/s')
    #read_draw_seasonal_4x3('q',1000,[0,1],'q','g/kg')
    #read_draw_seasonal_4x3('z','zonal',[-105,15],'ano_z','gpm')
    read_draw_seasonal_4x3('u',1,[-35,5],'U','km/s')
    #read_draw_seasonal_4x3('Q1',1,[-0.07,0.01],'Q1','W/kg')

def read_draw_seasonal_4x3(varname,scale,cnlev,label,unit):
    nrow = 4
    ncol = 3
    bmlo = 0.4

    f = xr.open_dataset('%s/ERA5_mon_%s_1979-2020.nc'%(path,varname))
    lat = f.latitude.data
    lon = f.longitude.data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]

    ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m

    if cnlev[0] < 0 :
        colr = cmaps.BlueDarkRed18(range(0,18,1))
        fcolors = colors.ListedColormap(np.vstack((colr[0:8,:],colr[10::,:])))
        print(fcolors.N)
    else:
        fcolors = cmaps.precip2_17lev
    cnlevels = np.arange(cnlev[0], cnlev[0]+cnlev[1]*(fcolors.N-1), cnlev[1])
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

    fig = plt.figure(figsize=(12,12),dpi=100)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree())) #sharex=True, sharey=True
    
    for nl in range(0,len(lev),1):
        ds = xr.open_dataset('%s/ERA5_mon_%s_1979-2020.nc'%(path,varname))
        da = ds[varname].sel(level=lev[nl],latitude=ilat,longitude=ilon,method="nearest").load()
        if scale=='zonal':
            #da = ds[varname].sel(level=lev[nl]).load()
            #var = da.groupby(da.time.dt.month).mean('time')/9.8
            #var.data = np.moveaxis(np.moveaxis(var.data,2,0)-var.data.mean(2), 0, 2)
            #var = var.sel(latitude=ilat,longitude=ilon,method="nearest").data
            var = da.groupby(da.time.dt.month).mean('time').data/9.8
            var = np.moveaxis(np.moveaxis(var,2,0)-var.mean(2), 0, 2)
        else:
            var = da.groupby(da.time.dt.month).mean('time').data*scale
        
        '''
        ds = xr.open_dataset('%s/ERA5_mon_z_1979-2020.nc'%path)
        da = ds['z'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
        var1 = da.groupby(da.time.dt.month).mean('time').data/9.8
        ds = xr.open_dataset('%s/ERA5_mon_u_1979-2020.nc'%path)
        da = ds['u'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
        uwnd = da.groupby(da.time.dt.month).mean('time').data
        ds = xr.open_dataset('%s/ERA5_mon_v_1979-2020.nc'%path)
        da = ds['v'].sel(level=lev[nl],longitude=ilon,latitude=ilat,method="nearest").load()
        vwnd = da.groupby(da.time.dt.month).mean('time').data
        var = dynamic_calc.calc_uv2vr_cfd(uwnd,vwnd,ilat,ilon) 
        '''
        del ds, da
        gc.collect()
        if lev[nl]==850:
            var = np.ma.array(var,mask=(
                np.broadcast_to(phis,var.shape)>1500))

        for nm in range(0,nrow,1):
            if nm == 0:
                shad = (var[0,:,:]+var[1,:,:]+var[11,:,:])/3.0
                #cont1 = (var1[0,:,:]+var1[1,:,:]+var1[11,:,:])/3.0
            else:
                shad = np.mean(var[(3*nm-1):(3*nm+2),:,:],axis=0)
                #cont1 = np.mean(var1[(3*nm-1):(3*nm+2),:,:],axis=0)
            axe = ax[nm][nl]
            axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            axe.set_title("(%s) %dhPa %s"%(numod[3*nm+nl],lev[nl],titls[nm]),
                fontsize=title_font,fontdict=font)

            print('min:%f ; max:%f'%(np.nanmin(shad),np.nanmax(shad)))
            shad = axe.contourf(ilon, ilat, shad, cnlevels,
                         transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
            
            #wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd1[::q_mis,::q_mis],vwnd1[::q_mis,::q_mis],
            #        pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="dimgray",
            #        width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

            #cont = axe.contour(ilon, ilat, cont1, np.arange(1000,15000,cnlvl2[nl]), 
            #             transform=ccrs.PlateCarree(), colors='darkviolet', linewidths=2)

            topo = axe.contour(ilon, ilat, phis, [1500,3000,4500],
                         transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)

            if nl == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.18, bmlo, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)
    #axe.quiverkey(wind, 0.92, bmlo-0.01, vcref[nl], r'$%d m/s$'%vcref[nl], labelpos='N',coordinates='figure')

    plt.figtext(0.02,bmlo-0.005, "%s (%s)"%(label,unit), fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir+"%s.png"%(varname), bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

