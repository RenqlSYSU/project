#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess, linecache, gc
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from renql import composite

title_font=14
label_font=10
lat_sp = 20
lon_sp = 30 #60 #
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

prefix = "fftadd"
radiu = 2.5
radius = 6
title= {'_%dlocal'%radius:'local',
        '_%doutside'%radius:'outside',
        '_%dtotal'%radius:'total',
        '':'All'}
suffixs = ["_%dlocal"%radius,"_%doutside"%radius]#'',"_%dtotal"%radius,
lev  = [850,500,250]
path = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev'
figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig/"
datapath = "/home/ys17-23/Extension2/renql/project/uor_track/mdata/"

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #

def main_run():
    #gene_grid()
    for ns in range(3):
        for nl in lev:
            calc_monthly_max_mean(suffixs[ns],nl,9,
                '%s/max_mean_10mwind_west%d%s.nc'%(datapath,nl,suffixs[ns]))
            
            if not os.path.exists('%s/max_mean_preci_west%d%s.nc'%(datapath,nl,suffixs[ns])):
                calc_monthly_max_mean(suffixs[ns],nl,13,
                    '%s/max_mean_preci_west%d%s.nc'%(datapath,nl,suffixs[ns]))
    '''
        ds = xr.open_dataset('%s/max_mean_preci_%d%s.nc'%(datapath,850,suffixs[ns]))
        ilat = ds['lat'].data
        ilon = ds['lon'].data
        dvars = ['maxv','numb','mean']
        for varname in ['preci']:#,'10mwind'
            if varname=='preci':
                cnlevs = [[50,900,50],[0,170,10],[10,350,20]]
                scale = 24
            else:
                cnlevs = [[2,27.5,1.5],[0,170,10],[0,17,1]]
                scale = 1
            for dvar,cnlev in zip(dvars,cnlevs):
                draw_seasonal_4x3(varname,suffixs[ns],dvar,cnlev,'%s %s'%(dvar,varname),
                    '%s/seasonal_%s_%s_%s.png'%(figdir,varname,dvar,suffixs[ns]),ilon,ilat,scale) 
    '''
def draw_seasonal_4x3(varname,suffix,dvar,cnlev,cblabel,figdir,ilon,ilat,scale):
    # varname = 'preci','10mwind', dvar='mean','maxv','numb'
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    month = ['DJF','MAM','JJA','SON']
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
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
        ds = xr.open_dataset('%s/max_mean_%s_%d%s.nc'%(
            datapath,varname,lev[nc],suffix))
        term = ds[dvar].data
        for nr in range(0,nrow,1):
            if nr == 0:
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
                if dvar=='maxv':
                    var = scale*np.max(np.array(
                        [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
                if dvar=='numb':
                    var = np.sum(np.array(
                        [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)/41.0
                if dvar=='mean':
                    var = scale*np.sum(np.array(
                        [term[0,:,:],term[1,:,:],term[11,:,:]]),axis=0)
                    term1 = ds['numb'].data
                    numb = np.sum(np.array(
                        [term1[0,:,:],term1[1,:,:],term1[11,:,:]]),axis=0)
                    var = np.where(numb>0,var/numb,0)
            else:
                uwnd1 = np.mean(uwnd[(3*nr-1):(3*nr+2),:,:],axis=0)
                if dvar=='maxv':
                    var = scale*np.max(term[(3*nr-1):(3*nr+2),:,:],axis=0)
                if dvar=='numb':
                    var = np.sum(term[(3*nr-1):(3*nr+2),:,:],axis=0)/41.0
                if dvar=='mean':
                    var = scale*np.sum(term[(3*nr-1):(3*nr+2),:,:],axis=0)
                    numb = np.sum(ds['numb'].data[(3*nr-1):(3*nr+2),:,:],axis=0)
                    var = np.where(numb>0,var/numb,0)
            print('%dhPa %s %s : min %f ; max %f'%(lev[nc],title[suffix],
                month[nr], var.min(), var.max()))
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%dhPa %s %s'%(lev[nc],title[suffix],month[nr]),
                fontsize=title_font, fontdict=font)

            cont = axe.contourf(ilon, ilat, var, cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
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

def calc_monthly_max_mean(suffix,nl,nint,fileout):
    # nint, intensity columns in fftadd file
    # nint=9,max10mwind; nint=10,average precip
    ds = xr.open_dataset("%s/grid_%.1f.nc"%(datapath,radiu))
    grid = ds['grid']
    dim = grid.shape
    print(dim)
    mean = np.zeros([12,dim[0],dim[1]],dtype=float)
    maxv = np.zeros([12,dim[0],dim[1]],dtype=float)
    numb = np.zeros([12,dim[0],dim[1]],dtype=float)

    ctime,clat,clon,cinte = composite.composite_time(
        "%s/%s_%d_1980-2020%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr,nint)
    nd = 0
    #for ct,clo,cla,cit in zip(ctime[0:100],clon[0:100],clat[0:100],cinte[0:100]):
    for ct,clo,cla,cit in zip(ctime,clon,clat,cinte):
        nd = nd + 1
        dist = np.square(grid.lat-cla)+np.square(grid.lon-clo)
        ind = np.unravel_index(np.argmin(dist.data), dist.shape)
        print('%d dist shape %s, ind %s'%(nd,str(dist.shape),str(ind)))
        numb[ct.month-1,ind[0],ind[1]] = numb[ct.month-1,ind[0],ind[1]] + 1
        mean[ct.month-1,ind[0],ind[1]] = mean[ct.month-1,ind[0],ind[1]] + cit
        if cit > maxv[ct.month-1,ind[0],ind[1]]:
            maxv[ct.month-1,ind[0],ind[1]] = cit
    
    print(numb)
    #mean = np.where(numb>0,mean/numb,0)
    ds = xr.Dataset(
            {
                "numb": (["month", "lat", "lon"], numb),
                "maxv": (["month", "lat", "lon"], maxv),
                "mean": (["month", "lat", "lon"], mean),
                },
            coords={
                "month": range(0,12,1), 
                "lat"  : (["lat"],grid.lat.data),
                "lon"  : (["lon"],grid.lon.data),
                },
            )
    ds.attrs["description"]='%d cyclone point'%(len(ctime))
    ds.to_netcdf(fileout,"w")

def gene_grid():
    lat = np.arange(lats,latn+0.1,radiu)
    lon = np.arange(lonl,lonr+0.1,radiu)
    var = np.zeros([len(lat),len(lon)],dtype=float)
    da = xr.DataArray(var, coords=[lat,lon], dims=["lat", "lon"])
    ds = da.to_dataset(name="grid")
    ds.to_netcdf("%s/grid_%.1f.nc"%(datapath,radiu),"w")

def monthly_contour(var,ilon,ilat,cnlev,figtitle,cblabel,figdir):
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
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
    nm = -1
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            nm = nm+1
            print('%s : min %f ; max %f'%(titls[nm],
                var[nm,:,:].min(), var[nm,:,:].max()))
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title(figtitle+" "+titls[nm],fontsize=title_font)

            cont = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            jets = axe.contour(ilon, ilat, uwnd[nm,:,:], [30,40,50], 
                 transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=1.5)
            
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.02,bmlo-0.005, cblabel,fontsize=title_font,
            horizontalalignment='left',verticalalignment='bottom')

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

