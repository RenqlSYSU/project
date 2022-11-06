#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess, linecache, gc
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from datetime import datetime

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

figdir = '/home/users/qd201969/uor_track/fig/'
datapath = "/work/scratch-pw2/renql/ERA5_hourly/wind10/ERA5_speed10_1hr_dec-jan"
years = range(1980,2021,1)

def main_run():
    lagd = 7
    lagh = 5
    lon=[60,62,58,110,90]
    lat=[50,45,35,45 ,35]
    #draw_daliy_speed4(i,lon,lat)
    #draw_hour_speed4(i,lon,lat)
    #draw_pdf_speed4(0,lagd,lagh,lon[0],lat[0])
    process_pool = Pool(processes=len(lat))
    results=[]
    #for nl in range(len(lat)):
    #    result=process_pool.apply_async(draw_daliy_speed4, 
    #            args=(nl,lon[nl],lat[nl],))
    #    results.append(result)
    for nl in range(len(lat)):
        result=process_pool.apply_async(draw_pdf_speed4, 
                args=(nl,lagd,lagh,lon[nl],lat[nl],))
        #result=process_pool.apply_async(draw_hour_speed4, 
        #        args=(nl+len(lat),lon[nl],lat[nl],))
        results.append(result)
    print(results) 
    print(results) 
    process_pool.close()
    process_pool.join() 
    print(results[0].get()) 

def draw_pdf_speed4(nk,lagd,lagh,lon,lat):
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(1,1,sharex=True)#, sharey=True
    for i,indx in enumerate([754,2940,5100,7350]):
        indx2=np.arange(indx-24*lagd,indx+24*(lagd+1),24)
        indx3=np.array([np.arange(d-lagh,d+lagh+1) 
            for d in indx2]).reshape((2*lagh+1)*(2*lagd+1))

        ds1  = xr.open_dataset("%s1980.nc"%(datapath))
        term = ds1['var1'][indx3].sel(lat=lat,lon=lon,method='nearest').data
        str_time = datetime.strptime(np.datetime_as_string(
                ds1['time'][indx].data,
                unit='h'), '%Y-%m-%dT%H').strftime('%b %d')
        for ny in range(1981,2021,1):
            ds1  = xr.open_dataset("%s%d.nc"%(datapath,ny))
            term = np.concatenate((term,
                ds1['var1'][indx3].sel(lat=lat,lon=lon,method='nearest').data))
            print("task[%d] %d time: %s"%(nk,ny,str_time),term.shape)
        
        var, xbin = np.histogram(term,20,density=True)
        ax.plot(xbin[1:],var,linewidth=2,label=str_time)
        
    ax.legend()
    ax.set_title("10m wind (%.2fN, %.2fE)"%(
            lat, lon),fontsize=title_font,fontdict=font)
    ax.set_xlabel('speed', fontsize=label_font, fontdict=font)
    fig.savefig("%s10wind_pdf_%dN_%dE_lagd%d_lagh%d.png"%(figdir,lat,lon,lagd,lagh),
            bbox_inches='tight',pad_inches=0.1)

def draw_hour_speed4(nk,lon,lat):
    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = fig.subplots(4,1,sharex=True)#, sharey=True
    for i,nh in enumerate([15,105,195,285]):
        ax = axe[i] 
        for ny in range(0,len(years)):
            ds = xr.open_dataset("%s%d.nc"%(datapath,years[ny]))
            stime = np.array([ds.time.dt.year.isin(years[ny]),
                    ds.time.dt.dayofyear.isin(nh)]).all(axis=0)
            var= ds['var1'].sel(time=stime,lat=lat,lon=lon,
                    method='nearest') 
            ax.plot(range(len(var)),var.data,linewidth=1)
            print("task[%02d]: %d time frame %d"%(nk,years[ny],len(var)))
        
        str_time = datetime.strptime(np.datetime_as_string(
            var.time[0]), '%Y-%m-%dT%H:00:00.000000000').strftime('%b %d')
        ax.set_title("10m wind %s (%.2fN, %.2fE)"%(str_time,
                var.lat, var.lon),fontsize=title_font,fontdict=font)
        #ax.set_xlabel('day', fontsize=label_font, fontdict=font)
        ax.set_ylabel('speed', fontsize=label_font, fontdict=font)
    fig.tight_layout()
    fig.savefig("%s10wind_hour_%dN_%dE.png"%(figdir,lat,lon),
            bbox_inches='tight',pad_inches=0.1)

def draw_daliy_speed4(nk,lon,lat):
    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = fig.subplots(4,1,sharex=True)#, sharey=True
    for i,nh in enumerate([0,6,12,18]):
        ax = axe[i] 
        for ny in range(0,len(years),4):
            ds = xr.open_dataset("%s%d.nc"%(datapath,years[ny]))
            stime = np.array([ds.time.dt.year.isin(years[ny]),
                    ds.time.dt.hour.isin(nh)]).all(axis=0)
            var= ds['var1'].sel(time=stime,lat=lat,lon=lon,
                    method='nearest') 
            ax.plot(range(len(var)),var.data,linewidth=1)
            print("task[%02d]: %d time frame %d"%(nk,years[ny],len(var)))

        ax.set_title("10m wind %02d:00 (%.2fN, %.2fE)"%(nh, 
                var.lat, var.lon),fontsize=title_font,fontdict=font)
        #ax.set_xlabel('day', fontsize=label_font, fontdict=font)
        ax.set_ylabel('speed', fontsize=label_font, fontdict=font)
    fig.tight_layout()
    fig.savefig("%s10wind_%d_%dN_%dE.png"%(figdir,nh,lat,lon),
            bbox_inches='tight',pad_inches=0.1)

def draw_daliy_speed(lon,lat):
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax  = fig.add_axes([0.05, 0.05, 0.9, 0.35])
    cmap_list = ['Greys','Purples','Greens','Reds']
    lgd = []
    
    for nh,clr in zip([0,6,12,18],cmap_list):
        cl = cm.get_cmap(clr)(range(30,156,3))

        ds = xr.open_dataset("%s1980.nc"%(datapath))
        stime = np.array([ds.time.dt.year.isin(1980),
                ds.time.dt.hour.isin(nh)]).all(axis=0)
        #term = np.argwhere(stime==True)
        #indx = np.array(term).reshape(len(term))
        var= ds['var1'].sel(time=stime,lat=lat,lon=lon,
                method='nearest') 
        print(var)
        
        line, = ax.plot(range(len(var)),var.data,linewidth=1,
                color=cl[0],label='%02d:00'%nh)
        lgd.append(line)
        for ny in range(1,len(years),4):
            ds = xr.open_dataset("%s%d.nc"%(datapath,years[ny]))
            stime = np.array([ds.time.dt.year.isin(years[ny]),
                    ds.time.dt.hour.isin(nh)]).all(axis=0)
            var= ds['var1'].sel(time=stime,lat=lat,lon=lon,
                    method='nearest') 
            ax.plot(range(len(var)),var.data,linewidth=2,color=cl[ny])
            print("%d time frame %d"%(years[ny],len(var)))

        ax.legend(handles=lgd,loc='upper right')
        ax.set_title("10m wind (%.2fN, %.2fE)"%(var.lat,var.lon),
                fontsize=title_font,fontdict=font)
        ax.set_xlabel('day',fontsize=label_font,fontdict=font)
        ax.set_ylabel('speed',fontsize=label_font,fontdict=font)
        fig.savefig("%s10wind_daily_%dN_%dE.png"%(figdir,lat,lon),
                bbox_inches='tight',pad_inches=0.1)

if __name__=='__main__':
    main_run()
