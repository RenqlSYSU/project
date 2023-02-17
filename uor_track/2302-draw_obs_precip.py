#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
import sys, os, subprocess
from datetime import datetime
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
path = '/home/ys17-23/Extension2/renql/'
lonl=10  #0  #
lonr=150#360#
lats=15 #20 #
latn=70 #90 #
ilat1 = np.arange(lats, latn+0.1, 1)
ilon1 = np.arange(lonl, lonr+0.1, 1)
lat_sp = 20
lon_sp = 30 #60 #
numod= [chr(i) for i in range(97,115)]

def main_run():
    outfile='%s/gpcp_monthly_preci.nc'%fileout
    if not os.path.exists(outfile):
        calc_gpcp_monthly(outfile)

    outfile='%s/gpm_monthly_preci.nc'%fileout
    if not os.path.exists(outfile):
        calc_gpm_monthly(outfile)
    ''' 
    ds = xr.open_dataset(outfile)
    ilon = ds.lon
    ilat = ds.lat
    draw_seasonal_preci_2x2(1,outfile,[0,17,1],'precip (mm/month)',
        '%s/preci.png'%(figdir),ilon,ilat)
    '''

    ds = xr.open_dataset('%s/clim_precip.nc'%fileout)
    lon = ds.longitude
    lat = ds.latitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    draw_seasonal_preci_4x3('precip (mm/month)',
        '%s/preci4x3.png'%(figdir),ilat,ilon)

def draw_seasonal_preci_4x3(cblabel,figdir,ilat,ilon):
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    month = ['DJF','MAM','JJA','SON']
    
    #cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    cnlevels = [0.1, 1, 3, 6, 10, 15, 25, 40, 60, 80, 100, 120, 150, 200, 250, 300, 350] #24h accumulated preci
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon1,lat=ilat1,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    uwndpath = '/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
    ds = xr.open_dataset(uwndpath)
    da = ds['u'].sel(level=200,longitude=ilon1,
        latitude=ilat1,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da
        
    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nc in range(3):
        if nc==0:
            ds = xr.open_dataset('%s/clim_precip.nc'%fileout)
            term = ds['tp'].sel(longitude=ilon,latitude=ilat).data
            title='ERA5'
        elif nc==1:
            ds = xr.open_dataset('%s/gpm_monthly_preci.nc'%fileout)
            term = ds['tp'].data.transpose(0,2,1)
            ilat = ds.lat
            ilon = ds.lon
            title='GPM'
        else:
            ds = xr.open_dataset('%s/gpcp_monthly_preci.nc'%fileout)
            term = ds['tp'].data
            ilat = ds.lat
            ilon = ds.lon
            title='GPCP'

        for nm in range(0,len(month),1):
            if nm == 0:
                term1 = np.mean(np.array(
                    [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                term1 = np.mean(term[(3*nm-1):(3*nm+2),:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
                
            print(term1.min())
            print(term1.max())
            axe = ax[nm][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title("(%s) %s %s"%(numod[nm*3+nc],title,month[nm]),
                fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, term1, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon1, ilat1, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            jets = axe.contour(ilon1, ilat1, uwnd1, [30,40,50], linestyles='solid', 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2.0)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.3, bmlo+0.005, 0.6, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.1,bmlo-0.005, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def calc_gpm_monthly(outfile):
    days = [31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,31   ,30   ,31   ]
    ds = xr.open_dataset('%s/GPM_mon/GPM_V06B_mon_200006-202105.nc'%path)
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    term = ds['precipitation'].sel(lat=ilat,lon=ilon)
    var = term.groupby(term.time.dt.month).mean('time')
    del term,ds

    for i in range(len(days)):
        var[i,:,:] = var[i,:,:]*days[i]*24
    var.attrs['units']='mm/month'
    ds2 = var.to_dataset(name='tp')
    ds2.to_netcdf(outfile,"w")

def calc_gpcp_monthly(outfile):
    days = [31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,31   ,30   ,31   ]
    ds = xr.open_dataset('%s/gpcp_precip.mon.mean.nc'%path)
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    term = ds['precip'].sel(lat=ilat,lon=ilon)
    var = term.groupby(term.time.dt.month).mean('time')
    del term,ds

    for i in range(len(days)):
        var[i,:,:] = var[i,:,:]*days[i]
    var.attrs['units']='mm/month'
    ds2 = var.to_dataset(name='tp')
    ds2.to_netcdf(outfile,"w")

def draw_seasonal_preci_2x2(scale,filname,cnlev,cblabel,figdir,ilon,ilat):
    nrow = 2 #6 #
    ncol = 2 #2 #
    bmlo = 0.4 #0.25 #
    month = ['DJF','MAM','JJA','SON']
    
    #cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    cnlevels = [0.1, 1, 3, 6, 10, 15, 25, 40, 60, 80, 100, 120, 150, 200, 250, 300, 350] #24h accumulated preci
    fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, 
        ncolors=fcolors.N,extend='both')
    
    ds = xr.open_dataset(filname)
    term = ds['tp'].sel(lon=ilon,lat=ilat).data.transpose(0,2,1)
    
    ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon1,lat=ilat1,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds

    uwndpath = '/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
    ds = xr.open_dataset(uwndpath)
    da = ds['u'].sel(level=200,longitude=ilon1,
        latitude=ilat1,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time').data
    del ds, da
        
    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0)))
    for nm in range(0,len(month),1):
        if nm == 0:
            term1 = scale*np.mean(np.array(
                [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
            uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
        else:
            term1 = scale*np.mean(term[(3*nm-1):(3*nm+2),:,:],axis=0)
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
        topo = axe.contour(ilon1, ilat1, phis, [1500,3000], 
             transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        jets = axe.contour(ilon1, ilat1, uwnd1, [30,40,50], linestyles='solid', 
             transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2.0)
        
        if nc == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nr == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.3, bmlo+0.07, 0.6, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.1,bmlo+0.06, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()
