#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys
from datetime import datetime
from renql import cyc_filter
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])
prefix = "fftadd"
suffix = sys.argv[1] 
radiu = int(sys.argv[2])
perc = int(sys.argv[3])
title= {'_local':'local',
        '_outside':'outside',
        '_total':'total',
        '':'All'}
behv = ["ALL" ,"PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]

fileout="/home/users/qd201969/uor_track/mdata/"
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
figdir = '/home/users/qd201969/uor_track/fig/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/precip/ERA5_precip_1hr_dec-jan" #1980.nc

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #
sy='clim'

def main_run():
    lag = 0

    #maxpreci_threshold(perc) 
    ds = xr.open_dataset("%smaxprecip1h_%dthreshold_month.nc"%(fileout,perc))
    ilon = ds.lon
    ilat = ds.lat
    thre = ds['threshold'].data
    print(thre)
    start_time = datetime.now()
    #draw_precip(thre*1000*24,[0,170,10],'%s'%sy,'max%dprecip (mm/day)'%perc
    #        ,figdir+'clim_precip',ilon,ilat,perc)
    #mask_precip(perc,var,0,850)
   
    process_pool = Pool(processes=3)
    results=[]
    for nl in range(len(lev)):
        result=process_pool.apply_async(mask_precip_amount, args=(perc,thre,lag,lev[nl],))
        results.append(result)
    print(results) 
    print(results[0].get()) 
    process_pool.close()
    process_pool.join() 
    print(results[0].get())
    print(start_time)
    print(datetime.now())
    print("%d precip lag %d: %s"%(perc,lag,title[suffix]))
    
    ''' 
    ds = xr.open_dataset("%s/clim_precip.nc"%(fileout))
    ilon = ds.longitude
    ilat = ds.latitude
    tp = ds['tp'].data
    #days = [31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,31   ,30   ,31   ]
    #for i in range(len(days)):
    #    tp[i,:,:] = tp[i,:,:]/days[i]
    #draw_precip(tp,[0,170,10],'%s'%sy,'precip (mm/month)'
    #        ,figdir+'clim_precip',ilon,ilat,perc)
    
    days = [31   ,28   ,31   ,30   ,31   ,30   ,31   ,31   ,30   ,31   ,30   ,31   ]
    ds = xr.open_dataset("%sclim_max%dprecip_event.nc"%(fileout,perc))
    ilon = ds.longitude
    ilat = ds.latitude
    var = ds['tp'].data
    var = tp
    #term1 = var
    #for i in range(len(days)):
    #    term1[i,:,:] = var[i,:,:]/days[i]
    #draw_precip(term1*30,[0,8.5,0.5],'%s extreme precip'%sy,'precip (h/30day)',figdir+'clim_precip',ilon,ilat,perc)

    for nl in lev:
        #ds = xr.open_dataset("%sclim_max%dprecip_%drad_lag%d_%d%s.nc"%(fileout,perc,radiu,lag,nl,suffix))
        ds = xr.open_dataset("%sclim_precip_%drad_lag%d_%d%s.nc"%(fileout,radiu,lag,nl,suffix))
        #ds = xr.open_dataset("%s%s_precip_%d%s.nc"%(fileout,sy,nl,suffix))
        #ds = xr.open_dataset("%s%s_precip_%d%s_%ddegree.nc"%(fileout,sy,nl,suffix,radiu))
        term = ds['tp'].data
        term = xr.where(var>0,(var-term)*100/var,0)
        draw_precip(term,[2,104,6],'%s %s %d'%(title[suffix],sy,nl),
                'precip percent','%sclim_precip%d%s'%(figdir,nl,suffix),ilon,ilat,perc)
        #ds = xr.open_dataset("%sclim_max%dprecip_6rad_lag%d_%d%s.nc"%(fileout,perc,lag,nl,suffix))
        #term = ds['tp'].data
        #term = xr.where(var>0,(var-term)*100/var,0)
        #ds = xr.open_dataset("%sclim_max%dprecip_7rad_lag%d_%d%s.nc"%(fileout,perc,lag,nl,suffix))
        #term2 = ds['tp'].data
        #term2 = xr.where(var>0,(var-term2)*100/var,0)
        #draw_precip(term2-term,[0,51,3],'7rad-6rad %s %s %d'%(title[suffix],sy,nl),
        #        'precip percent','%sclim_precip_diff%d%s'%(figdir,nl,suffix),ilon,ilat,perc)
    '''
def maxpreci_threshold(perc):
    ds  = xr.open_dataset(datapath+"1980.nc")
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    thre = np.empty( [12,len(ilat),len(ilon)],dtype=float ) 
    
    for nm in range(0,12,1):
        var = ds['tp'].sel(time=np.array(
            [ds.time.dt.month.isin(nm+1),ds.time.dt.year.isin(1980)]
            ).all(axis=0),longitude=ilon,latitude=ilat).data
        for ny in range(1981,2021,1):
            ds1  = xr.open_dataset("%s%d.nc"%(datapath,ny))
            term = ds1['tp'].sel(time=np.array(
                [ds1.time.dt.month.isin(nm+1),ds1.time.dt.year.isin(ny)]
                ).all(axis=0),longitude=ilon,latitude=ilat).data
            var = np.concatenate((var, term))
            print("month %2d %d: "%(nm,ny), var.shape)
        thre[nm,:,:] = np.percentile(var,perc,axis=0)
    
    da = xr.DataArray(thre, coords=[range(0,12,1),ilat,ilon], 
            dims=["month","lat","lon"])
    ds2 = da.to_dataset(name='threshold')
    ds2.to_netcdf("%smaxprecip1h_%dthreshold_month.nc"%(fileout,perc),"w")

def mask_precip_amount(perc,thre,lag,nl):
    ds  = xr.open_dataset(datapath+"1980.nc")
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var  = np.empty( [12,len(ilat),len(ilon)],dtype=float )  
    
    ctime,clat,clon = composite_time("%s%s_%d_1980-2020%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr)
    for ny in range(1980,2021,1):
        print('task [%d]:%d'%(nl,ny))
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,longitude=ilon,latitude=ilat)
        
        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            indx = np.argwhere(term.time.data==ct.data)[0][0]
            term[indx-lag:(indx+1+lag),:,:] = term[indx-lag:(indx+1+lag),:,:].where(
                (np.square(term.longitude-clo)+
                np.square(term.latitude-cla))>(radiu*radiu), 0)
        var = var + term.groupby(term.time.dt.month).sum('time') 
    var = var*1000/41
    var.attrs['units']='mm/month'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%sclim_precip_%drad_lag%d_%d%s.nc"%(fileout,radiu,lag,nl,suffix),"w")

def mask_precip_extreme(perc,thre,lag,nl):
    ds  = xr.open_dataset(datapath+"1980.nc")
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    var  = np.empty( [12,len(ilat),len(ilon)],dtype=float )  
    '''
    for ny in range(1980,2021,1):
        print(ny)
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,longitude=ilon,latitude=ilat)
        #for nm in range(12):
        #    term.loc[dict(time=term.time.dt.month.isin(nm+1))] = xr.where(
        #        term.sel(time=term.time.dt.month.isin(nm+1))>thre[nm,:,:],1,0)
        var = var + term.groupby(term.time.dt.month).sum('time') 
    #var = var*1000/41
    #var.attrs['units']='mm/month'
    var = var/41
    var.attrs['units']='event'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%sclim_max%dprecip_event.nc"%(fileout,perc),"w")
    '''
    ctime,clat,clon = composite_time("%s%s_%d_1980-2020%s"%(
        path,prefix,nl,suffix),lats,latn,lonl,lonr)
    for ny in range(1980,2021,1):
        print('task [%d]:%d'%(nl,ny))
        ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
        term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
                ,longitude=ilon,latitude=ilat)
        #term = xr.where(term>=perc,1,0)
        for nm in range(12):
            term.loc[dict(time=term.time.dt.month.isin(nm+1))] = xr.where(
                term.sel(time=term.time.dt.month.isin(nm+1))>thre[nm,:,:],1,0)
        
        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            indx = np.argwhere(term.time.data==ct.data)[0][0]
            term[indx-lag:(indx+1+lag),:,:] = term[indx-lag:(indx+1+lag),:,:].where(
                (np.square(term.longitude-clo)+
                np.square(term.latitude-cla))>(radiu*radiu), 0)
        #term = term/term.time.dt.days_in_month
        var = var + term.groupby(term.time.dt.month).sum('time') 
    #var = var*1000/41
    #var.attrs['units']='mm/month'
    var = var/41
    var.attrs['units']='event'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%sclim_max%dprecip_%drad_lag%d_%d%s.nc"%(fileout,perc,radiu,lag,nl,suffix),"w")
    #ds1.to_netcdf("%sclim_precip_%drad_lag%d_%d%s.nc"%(fileout,radiu,lag,nl,suffix),"w")

def mask_precip_oneyear(ny):
    ds  = xr.open_dataset("%s%d.nc"%(datapath,ny))
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    term = ds['tp'].sel(time=ds.time.dt.year.isin(ny)
            ,longitude=ilon,latitude=ilat)
    var = term.groupby(term.time.dt.month).sum('time') 
    var.attrs['units']='mm/month'
    ds1 = var.to_dataset(name='tp')
    ds1.to_netcdf("%s%d_precip.nc"%(fileout,ny),"w")

    for nl in lev:
        ctime,clat,clon = composite_time("%s%s_%d_1980-2020%s"%(
            path,prefix,nl,suffix),lats,latn,lonl,lonr)

        ctime1 = ctime.where(ctime.dt.year.isin(ny),drop=True)
        clat1 = clat.where(ctime.dt.year.isin(ny),drop=True)
        clon1 = clon.where(ctime.dt.year.isin(ny),drop=True)
        print(ctime1)
        for ct,clo,cla in zip(ctime1,clon1,clat1):
            print(ct)
            term.loc[ct,:,:] = term.sel(time=ct).where(
                (np.square(term.longitude-clo)+np.square(term.latitude-cla))>25, 0)
        var.data = term.groupby(term.time.dt.month).sum('time').data
        ds1 = var.to_dataset(name='tp')
        ds1.to_netcdf("%s%d_precip_%d%s.nc"%(fileout,ny,nl,suffix),"w")

def draw_precip(var,cnlev,figtitle,cblabel,figdir,ilon,ilat,perc):
    sample = np.array([30504,27816,30504,29520,30504,29520,30504,30504,29520,30504,29520,30504])
    sample = sample*(1-perc)/41
    print(sample)
    
    lat_sp = 20
    lon_sp = 30 #60 #
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.37 #0.25 #
    title_font=14
    label_font=10
    titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    
    ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    
    #cnlevels = np.arange(cnlev[0], cnlev[1], cnlev[2])
    cnlevels = [0.1, 1, 3, 6, 10, 15, 25, 40, 60, 80, 100, 120, 150, 200, 250, 300, 350] #24h accumulated preci
    #cnlevels = [0.1, 0.5, 1, 2, 3, 5, 8, 12, 16, 20, 25, 30, 40, 50, 70, 100, 150] #houly preci
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
            if perc<1 :
                var[nm,:,:] = (sample[nm]-var[nm,:,:])*100/sample[nm]
                print('percent')
            print(var[nm,:,:].min())
            print(var[nm,:,:].max())
            axe = ax[nr][nc]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title(figtitle+" "+titls[nm],fontsize=title_font)

            cont = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                 transform=ccrs.PlateCarree(),cmap=fcolors,
                 extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                 linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd[nm,:,:], [30,40,50], 
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
        position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
        plt.figtext(0.02,bmlo-0.005, cblabel,fontsize=title_font,
                horizontalalignment='left',verticalalignment='bottom')
    
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir, bbox_inches='tight',pad_inches=0.01)

def composite_time(filname,flats,flatn,flonl,flonr):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    ctime=[]
    clat =[]
    clon =[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            for nl in range(0,num,1):
                line = ff.readline()
                if prefix in ['ffadd','fftadd']:
                    data = list(map(float,line.strip().replace(" &","").split(" ")))
                else:
                    data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    ctime.append(datetime.strptime(str(int(data[0])),'%Y%m%d%H'))
                    clat.append(data[2])
                    clon.append(data[1])

        line = ff.readline()
    ff.close()
    ctime = xr.DataArray(ctime)
    clat = xr.DataArray(clat)
    clon = xr.DataArray(clon)
    return ctime,clat,clon

if __name__=='__main__':
    main_run()

