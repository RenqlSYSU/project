#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess
from datetime import datetime
from renql import cyc_filter
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig/"
fileout = "/home/ys17-23/Extension2/renql/project/uor_track/mdata/"
path = '/home/ys17-23/Extension2/renql/ERA5_mon'
path2= '/home/ys17-23/Extension2/renql/ERA5-1HR-lev/statistic'
lev  = [850,500,250]
title= {'_6local':'local',
        '_6outside':'outside',
        '_total':'total',
        '':'All'}
#suffix = '_6local'#,'_total']'_6outside'#,''#
suffix = '_6outside'
#suffix = ''
if suffix in ['_6local','local','remote']:
    lonl=50  #0  #
    lonr=150#360#
    lats=15 #20 #
    latn=55 #90 #
    bmlo = 0.37 #0.25 #
if suffix in ['_6outside','outside']:
    lonl=20  #0  #
    lonr=140#360#
    lats=15 #20 #
    latn=60 #90 #
    bmlo = 0.37 #0.25 #
cnlev2 = np.hstack((np.arange(1,8.5,2.5),np.arange(13.5,50,5)))
dash = [8.5,60]
if suffix in ['']:
    lonl=0  #0  #
    lonr=150#360#
    lats=15 #20 #
    latn=70 #90 #
    bmlo = 0.37 #0.25 #

#lonl=0  #0  #
#lonr=150#360#
#lats=15  #
#latn=70 #
lat_sp = 20
lon_sp = 30 #60 #
# grid used for contribution 
#ilat = np.arange(latn, lats-0.1,-0.25)
#ilon = np.arange(lonl, lonr+0.1, 0.25)
# grid used for max preci induced by cyclone
ilat = np.arange(lats, latn+0.1, 2.5)
ilon = np.arange(lonl, lonr+0.1, 2.5)
radiu = 6
perc = 99
numod= [chr(i) for i in range(97,115)]

ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds

def main_run():
    ds = xr.open_dataset("%s/clim_precip.nc"%(fileout))
    var = ds['tp'].sel(longitude=ilon,latitude=ilat).data
    print(var)
    draw_seasonal_contour4x3_diff(var,'tp','%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu),
        suffix,[2,104,6],'maxpreci (%)',
        '%s/max%dprecip%s_contri_4x3.png'%(figdir,perc,suffix),ilon,ilat)
    
    '''
    # draw extreme contribution 
    ds = xr.open_dataset("%s/clim_max%dprecip_event.nc"%(fileout,perc))
    var = ds['tp'].sel(longitude=ilon,latitude=ilat).data
    print(var)
    draw_seasonal_contour4x3_diff(var,'tp','%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu),
        suffix,[2,104,6],'maxpreci (%)',
        '%s/max%dprecip%s_contri_4x3.png'%(figdir,perc,suffix),ilon,ilat)
    draw_seasonal_contour4x3_diff(var,'event','%s/clim_%.1fmax10mwind_%drad_lag0'%(fileout,perc,radiu),
        suffix,[2,104,6],'max10mwind (%)',
        '%s/%.1fmax10mwind%s_contri_4x3.png'%(figdir,perc,suffix),ilon,ilat)
    # draw max weather induced by cyclone
    draw_seasonal_4x3(1,'%s/max_mean_preci'%fileout,suffix,
        [2,53,3],'maxprecip (mm/h)','%s/max_precip%s.png'%(figdir,suffix),ilon,ilat)
    draw_seasonal_4x3(1,'%s/max_mean_10mwind'%fileout,suffix,
        [2,27.5,1.5],'max10mwind (m/s)','%s/max_10mwind%s.png'%(figdir,suffix),ilon,ilat)

    # threshold
    draw_seasonal_2x2(1,'%s/max10mwind_99.0threshold_month.nc'%fileout,
        [2,19,1],'max10mwind (m/s)','%s/threshold_10mwind.png'%(figdir),ilon,ilat)
    draw_seasonal_2x2(1000,'%s/maxprecip1h_99threshold_month.nc'%fileout,
        [0,8.5,0.5],'maxprecip (mm/h)','%s/threshold_precip.png'%(figdir),ilon,ilat)
    '''

def draw_seasonal_contour4x3_diff(var,varname,filname,suffix,cnlev,cblabel,figdir,ilon,ilat):
    nrow = 4 #6 #
    ncol = 3 #2 #
    month = ['DJF','MAM','JJA','SON']
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nc in range(0,ncol,1):
        ds = xr.open_dataset("%s_%d%s.nc"%(filname,lev[nc],suffix))
        if varname=='tp':
            term = ds[varname].sel(longitude=ilon,latitude=ilat).data
        else:
            term = ds[varname].sel(lon=ilon,lat=ilat).data
        files = '%s/ff_%d_1980-2020%s_stat.nc'%(path2,lev[nc],suffix)
        ds = xr.open_dataset(files) 
        lat = ds.lat
        lon = ds.long
        ilon1 = lon[(lon>=lonl) & (lon<=lonr)]
        ilat1 = lat[(lat>=lats) & (lat<=latn)]
        print(files)
        uwnd = read_stat(ds,'tden',lev[nc],ilat1,ilon1)
        
        for nr in range(0,nrow,1):
            if nr == 0:
                term1 = term[0,:,:]+term[1,:,:]+term[11,:,:]
                var1 = var[0,:,:]+var[1,:,:]+var[11,:,:]
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                print(term[(3*nr-1):(3*nr+2),:,:].shape)
                term1 = np.sum(term[(3*nr-1):(3*nr+2),:,:],axis=0)
                var1 = np.sum(var[(3*nr-1):(3*nr+2),:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nr-1):(3*nr+2),:,:],axis=0)
            term1 = xr.where(var1>0,(var1-term1)*100/var1,0)
            if lev[nc]==850:
                term1 = np.ma.array(term1, mask=(phis>1500))
            print(term1.min())
            print(term1.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title("(%s) %dhPa %s %s"%(numod[3*nr+nc],lev[nc],month[nr],
                title[suffix]),fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, term1, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            #jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
            #     transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            #line2 = axe.contour(ilon1, ilat1, uwnd1, cnlev2, linestyles='solid', 
            #     transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            #line4 = axe.contour(ilon1, ilat1, uwnd1, dash, linestyles='dashed', 
            #     transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.3, bmlo+0.005, 0.6, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.1,bmlo-0.005, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def draw_seasonal_4x3(scale,filname,suffix,cnlev,cblabel,figdir,ilon,ilat):
    nrow = 4 #6 #
    ncol = 3 #2 #
    month = ['DJF','MAM','JJA','SON']
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nc in range(0,ncol,1):
        ds = xr.open_dataset('%s_%d%s.nc'%(filname,lev[nc],suffix))
        term = ds['maxv'].sel(lon=ilon,lat=ilat).data
        
        files = '%s/ff_%d_1980-2020%s_stat.nc'%(path2,lev[nc],suffix)
        ds = xr.open_dataset(files) 
        lat = ds.lat
        lon = ds.long
        ilon1 = lon[(lon>=lonl) & (lon<=lonr)]
        ilat1 = lat[(lat>=lats) & (lat<=latn)]
        print(files)
        uwnd = read_stat(ds,'tden',lev[nc],ilat1,ilon1)
        
        for nr in range(0,nrow,1):
            if nr == 0:
                term1 = scale*np.max(np.array(
                    [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                term1 = scale*np.max(term[(3*nr-1):(3*nr+2),:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nr-1):(3*nr+2),:,:],axis=0)
            
            if lev[nc]==850:
                term1 = np.ma.array(term1, mask=(phis>1500))
            
            print(term1.min())
            print(term1.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title("(%s) %dhPa %s %s"%(numod[3*nr+nc],lev[nc],month[nr],
                title[suffix]),fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, term1, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            line2 = axe.contour(ilon1, ilat1, uwnd1, cnlev2, linestyles='solid', 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            line4 = axe.contour(ilon1, ilat1, uwnd1, dash, linestyles='dashed', 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.3, bmlo+0.005, 0.6, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.1,bmlo-0.005, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def read_stat(f,varname,nl,ilat1,ilon1):
    var = f[varname].sel(long=ilon1,lat=ilat1).load()
    if varname == 'mten':
        var.data = var.data*24
    if varname[-3:] != 'den':
        tden = f['tden'].sel(long=ilon1,lat=ilat1).data
        mask = tden < 1.0
        var.data = np.ma.array(var.data,mask=mask)
        del mask, tden
    
    if nl==850:
        ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
        phis1 = ds['PHIS'].sel(lon=ilon1,lat=ilat1,method="nearest").load()
        phis1 = phis1/9.8 # transfer from m2/s2 to m
        del ds
        var.data = np.ma.array(var.data,mask=(
            np.broadcast_to(phis1, var.shape)>1500))
    
    print('%s: min:%f ; max:%f'%(var.long_name,
        np.nanmin(var.data),np.nanmax(var.data)))
    return var.data
    
def draw_seasonal_2x2(scale,filname,cnlev,cblabel,figdir,ilon,ilat):
    nrow = 2 #6 #
    ncol = 2 #2 #
    bmlo = 0.4 #0.25 #
    month = ['DJF','MAM','JJA','SON']
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    ds = xr.open_dataset(filname)
    term = ds['threshold'].sel(lon=ilon,lat=ilat).data
    
    uwndpath = '/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
    ds = xr.open_dataset(uwndpath)
    da = ds['u'].sel(level=200,longitude=ilon,
        latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da
        
    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nm in range(0,len(month),1):
        if nm == 0:
            term1 = scale*np.max(np.array(
                [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
            uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
        else:
            term1 = scale*np.max(term[(3*nm-1):(3*nm+2),:,:],axis=0)
            uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
            
        print(term1.min())
        print(term1.max())
        nr = int(nm/2)
        nc = nm-2*nr
        axe = ax[nr][nc]
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                , linewidth=0.8, zorder=1)
        axe.set_title("(%s) %s"%(numod[nm],month[nm]),
            fontsize=title_font,fontdict=font)

        cont = axe.contourf(ilon, ilat, term1, cnlevels, 
             transform=ccrs.PlateCarree(),cmap=fcolors,
             extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000], 
             transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], linestyles='solid', 
             transform=ccrs.PlateCarree(),colors='r',linewidths=2.0)
        
        if nc == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nr == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.3, bmlo+0.08, 0.6, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.1,bmlo+0.075, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()
