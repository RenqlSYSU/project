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
lev  = [850,500,250]
dbox = False
title= {'_local':'local',
        '_outside':'outside',
        '_total':'total',
        '':'All'}
suffixs = ['_local','_outside']#,'_total',''
lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #
lat_sp = 20
lon_sp = 30 #60 #
flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 110 #int(sys.argv[5])
radiu = 6
perc = 99

def main_run():
    #calc_associated_weather()
    ''' 
    ds = xr.open_dataset("%s/clim_precip.nc"%(fileout))
    ilon = ds.longitude
    ilat = ds.latitude
    var = ds['tp'].data
    draw_annual_contour3x3(var,'tp','%s/clim_precip_%drad_lag0'%(fileout,radiu),
        [2,104,6],'precip (%)',
        '%s/precip_contribution_3x3.jpg'%(figdir),ilon,ilat,suffixs[0:3])
    for ns in range(len(suffixs)):
        draw_seasonal_contour4x3(var,'tp','%s/clim_precip_%drad_lag0'%(fileout,radiu),
            suffixs[ns],[2,104,6],'preci (%)',
            '%s/precip%s_contri_4x3.jpg'%(figdir,suffixs[ns]),ilon,ilat)
    ''' 
    
    ds = xr.open_dataset("%s/clim_max%dprecip_event.nc"%(fileout,perc))
    var = ds['tp'].data
    #draw_annual_contour3x3(var,'tp','%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu),
    #    [2,104,6],'max precip (%)',
    #    '%s/max%dprecip_contribution_3x3.jpg'%(figdir,perc),ilon,ilat,suffixs[0:3])
    for ns in range(len(suffixs)):
        ilat = np.arange(latn, lats-0.1,-0.25)
        ilon = np.arange(lonl, lonr+0.1, 0.25)
        draw_seasonal_contour4x3(var,'tp','%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu),
            suffixs[ns],[2,104,6],'maxpreci (%)',
            '%s/max%dprecip%s_contri_4x3.jpg'%(figdir,perc,suffixs[ns]),ilon,ilat)
    
    ds = xr.open_dataset("%sclim_%.1fmax10mwind_event.nc"%(fileout,perc))
    var = ds['event'].data
    #draw_annual_contour3x3(var,'event','%s/clim_%.1fmax10mwind_%drad_lag0'%(fileout,perc,radiu),
    #    [2,104,6],'max10mwind (%)',
    #    '%s/%.1fmax10mwind_contribution_3x3.jpg'%(figdir,perc),ilon,ilat,suffixs[0:3])
    for ns in range(len(suffixs)):
        ilat = np.arange(latn, lats-0.1,-0.25)
        ilon = np.arange(lonl, lonr+0.1, 0.25)
        draw_seasonal_contour4x3(var,'event','%s/clim_%.1fmax10mwind_%drad_lag0'%(fileout,perc,radiu),
            suffixs[ns],[2,104,6],'max10mwind (%)',
            '%s/%.1fmax10mwind%s_contri_4x3.jpg'%(figdir,perc,suffixs[ns]),ilon,ilat)
    #draw_annual_contour3x2([2,104,6],'percent',
    #    '%s/max%dprecip_10mwind_contri_3x2.jpg'%(figdir,perc),ilon,ilat)

def calc_associated_weather():
    for suffix in suffixs[0:2]: 
        com = "python ~/uor_track/2203-calc_maximum10mwind2.py %s %d %d"\
                %(suffix,radiu,perc)
        print(com)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()
        
        com = "python ~/uor_track/2203-calc_clim_precip_mpool.py %s %d %d"\
                %(suffix,radiu,perc)
        print(com)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()

def draw_annual_contour3x2(cnlev,cblabel,figdir,ilon,ilat):
    nrow = 3 #6 #
    ncol = 2 #2 #
    bmlo = 0.32 #0.25 #
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.mean('time').data
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nc in range(0,ncol,1):
        if nc == 0:
            title1 = 'Extreme Preci'
            filname = '%s/clim_max%dprecip_%drad_lag0'%(fileout,perc,radiu)
            ds = xr.open_dataset("%sclim_max%dprecip_event.nc"%(fileout,perc))
            var = ds['tp'].data
            varname = 'tp'
        if nc == 1:
            title1 = 'Extreme 10mwind'
            filname = '%s/clim_%.1fmax10mwind_%drad_lag0'%(fileout,perc,radiu)
            ds = xr.open_dataset("%sclim_%.1fmax10mwind_event.nc"%(fileout,perc))
            var = ds['event'].data
            varname = 'event'
        var = np.sum(var,axis=0) 
        for nr in range(0,nrow,1):
            ds = xr.open_dataset("%s_%d%s.nc"%(filname,lev[nr],suffixs[2]))
            term = ds[varname].data
            term = np.sum(term,axis=0) 
            term = xr.where(var>0,(var-term)*100/var,0)
            #term = np.mean(term,axis=0)
            print(term.min())
            print(term.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%s %dhPa'%(title1, lev[nr]),
                fontsize=title_font, fontdict=font)

            cont = axe.contourf(ilon, ilat, term, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            if dbox:
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                     linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd, [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.99, bmlo+0.1, 0.01, 0.5]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
    cb.set_label(label=cblabel, size=title_font) #, weight='bold'
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def draw_seasonal_contour4x3(var,varname,filname,suffix,cnlev,cblabel,figdir,ilon,ilat):
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    month = ['DJF','MAM','JJA','SON']
    
    ds = xr.open_dataset('/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da
    
    ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nc in range(0,ncol,1):
        ds = xr.open_dataset("%s_%d%s.nc"%(filname,lev[nc],suffix))
        if varname=='tp':
            term = ds[varname].sel(longitude=ilon,latitude=ilat).data
        else:
            term = ds[varname].sel(lon=ilon,lat=ilat).data
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
            print(term1.min())
            print(term1.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%dhPa %s %s'%(lev[nc],title[suffix],month[nr]),
                fontsize=title_font, fontdict=font)

            cont = axe.contourf(ilon, ilat, term1, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            if dbox:
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                     linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.02,bmlo+0.005, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def draw_annual_contour3x3(var,varname,filname,cnlev,cblabel,figdir,ilon,ilat,suffix):
    nrow = 3 #6 #
    ncol = 3 #2 #
    bmlo = 0.45 #0.25 #
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.mean('time').data
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    var = np.sum(var,axis=0) 
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            ds = xr.open_dataset("%s_%d%s.nc"%(filname,lev[nr],suffix[nc]))
            term = ds[varname].data
            term = np.sum(term,axis=0) 
            term = xr.where(var>0,(var-term)*100/var,0)
            print(term.min())
            print(term.max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%dhPa %s'%(lev[nr],title[suffix[nc]]),
                fontsize=title_font, fontdict=font)

            cont = axe.contourf(ilon, ilat, term, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            if dbox:
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                     linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd, [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    if (nrow/ncol) >= 2.0: 
        position = fig.add_axes([0.99, bmlo+0.05, 0.01, 0.6]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
        cb.set_label(label=cblabel, size=title_font) #, weight='bold'
    else:
        position = fig.add_axes([0.2, bmlo+0.04, 0.7, 0.01]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
        plt.figtext(0.02,bmlo+0.03, cblabel,fontsize=title_font,
                horizontalalignment='left',verticalalignment='bottom')
    
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()
