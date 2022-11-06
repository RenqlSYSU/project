#!/usr/bin/env python
'''
1. when cyclone reach the defined region, read datetime
    , one cyclone only one time be read
2. calc the cyclone number time series ([len(lev),len(behv),len(months),len(years)] 
3. plot the time series and store the data

2021-11-08
renql
'''

import xarray as xr
import numpy as np
import pandas as pd
import sys, os, subprocess, linecache, gc
from datetime import datetime
from scipy import stats
from renql import monthly_calc, life_intensity, cyc_filter
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import colors

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

def draw_3var_distri():
    #varname = ["close1 (%)","close2 (%)","GeopoHeight (m)"]
    varname = ["10mWind (m/s)","MaxRain (mm/d)","RainTime (%)"]
    nrow = len(lev)
    ncol = 3
    bmlo = 0.4
    #xbin = [np.arange(0,101,5), # num=(end-start)/inv
    #        np.arange(0,101,5), 
    #        np.arange(0,8001,400)]
    xbin = [np.arange(2,22.1,1), # num=(end-start)/inv
            np.arange(0,50.1,2.5), 
            np.arange(0,101,5)]
    var = np.empty( [len(varname),len(behv),len(xbin[0])-1],dtype=float )   
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for nl in range(0,len(lev),1):
        #if nl==0 : # for geopotential
        #    xbin[2]=np.arange(1200,1701,25)
        #elif nl==1 :
        #    xbin[2]=np.arange(5000,6001,50)
        #elif nl==2 :
        #    xbin[2]=np.arange(9500,11001,75)

        numb = []
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"
        for nb in range(1,len(behv),1):
            if nb == 0:
                filname1 = filname
            else:
                suffix = str(option[nb])+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
                filname1 = filname+"_"+suffix
            
            life1, inte1, dist1, numb1 = life_intensity.calc_close_cyclone(
                    filname1,flats=20,flatn=50,flonl=105,flonr=130)
            print("%dhPa %s :%d"%(lev[nl],behv[nb],numb1))
            var[0,nb,:],term = np.histogram(life1,xbin[0])
            var[1,nb,:],term = np.histogram(inte1,xbin[1])
            var[2,nb,:],term = np.histogram(dist1,xbin[2])
            var[:,nb,:] = var[:,nb,:]*100/numb1
            numb.append(numb1)
                
        for nv in range(0,3,1):
            ax[nl][nv].grid(True, which="both", color='grey', linestyle='--', linewidth=1)
            ax[nl][nv].yaxis.set_major_formatter(mtick.FormatStrFormatter('%i'))
            for nb in range(1,len(behv),1):
        #        print(var[nv,nb,:])
                ax[nl][nv].plot(xbin[nv][1:],var[nv,nb,:],linewidth=2)
            
            if nl==2 :
                ax[nl][nv].set_xlabel(varname[nv],fontsize=label_font-4,fontdict=font)
            if nv==0 :
                ax[nl][nv].set_ylabel("%dhPa percent"%lev[nl],fontsize=label_font-4,fontdict=font)

    ax[0][2].legend(behv[1:len(behv)], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sclse_pinte_%dh_total"%(figdir,time))

def draw_ts(var,figdir):
    plinsty = ["solid",(0,(4,2)),(0,(1,2))] # change with lev
    pcolor  = ["k","r","b"] # change with option 

    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = plt.axes()
    axe.set_title("%d-%dE, %d-%dN"%(flonl,flonr,flats,flatn),fontsize=title_font,fontdict=font)
    axe.set_xlim(1, 12)
    axe.set_ylim(0, np.max(var))
    axe.set_xlabel("month",fontsize=label_font,fontdict=font)
    axe.set_ylabel("number of tracks",fontsize=label_font,fontdict=font)
    axe.tick_params(axis='both', which='major', labelsize=label_font)
    for nl in range(len(lev)):
        for nop in range(len(option)):
            axe.plot(months,var[nl,nop,:],linewidth=1.2+nl/2.0,color=pcolor[nop],linestyle=plinsty[nl],
                label=str(lev[nl])+" "+behv[nop]) #
    axe.legend(ncol=3,fontsize=label_font)
    fig.savefig(figdir+"annual_ts.png", bbox_inches='tight',pad_inches=0.01)

if len(sys.argv) < 2 :
    flats = 25 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 110 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "ff"
else:
    flats = int(sys.argv[1])
    flatn = int(sys.argv[2])
    flonl = int(sys.argv[3])
    flonr = int(sys.argv[4])
    time = int(sys.argv[5])

option = [2,0,5]
behv   = ["total","local","outside"]
figdir ="/home/users/qd201969/uor_track/fig/"
path = '/home/users/qd201969/ERA5-1HR-lev/'
lev  = [850,500,250]
months = range(1,13,1)
timefilt = 1
draw_hist = 0
draw_annual_ts = 0

if timefilt == 1:
    for nop in range(0,len(option)):
        suffix = str(option[nop])+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
        com = "sh ~/uor_track/control_era5_1hr_track.sh %s 2 1 _%s"\
            %(prefix,suffix)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()
        #for nl in range(len(lev)):
        #    filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix
            #numb = cyc_filter.time_filt(filname,time)
        # box filter 
        #com = "sh ~/uor_track/control_era5_1hr_track.sh %s 3 1 %d %d %d %d %d"\
        #    %(prefix,option[nop],flats,flatn,flonl,flonr)
        # calc statistics and draw figure 
#draw_3var_distri()

if draw_annual_ts == 1:
    var = np.zeros( [len(lev),3,12], dtype=float )
    for nl in range(len(lev)):
        for nop in range(len(option)):
            suffix = str(option[nop])+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
            filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix
            var[nl,nop,:] = monthly_calc.calc_month(filname,time)/41.0
    draw_ts(var,figdir)

if draw_hist == 1:
    nrow = len(lev)
    ncol = 2 
    bmlo = 0.1
    title_font=18
    label_font=14
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol, sharex=True, sharey=True)
    for nl in range(len(lev)):
        for nop in range(1,len(option)):
            suffix = str(option[nop])+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
            filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix
            patches = life_intensity.hist_life_intensity(filname, \
                    ax=ax[nl,nop-1], title="%d %s"%(lev[nl],behv[nop]),\
                    flats=flats,flatn=flatn,flonl=flonl,flonr=flonr)
            if nl == 2:
                ax[nl,nop-1].set_xlabel("lifetime (days)",fontsize=label_font,fontdict=font)

    ax[1,0].set_ylabel("Mean intensity ($10^{-5} s^{-1}$)",fontsize=label_font,fontdict=font)
    position = fig.add_axes([0.96, bmlo+0.05, 0.02, 0.8]) #left, bottom, width, height
    cb = plt.colorbar(patches, cax=position ,orientation='vertical')#, shrink=.9)
    fig.tight_layout(rect=(0,bmlo,0.94,1)) #w_pad=0.5,h_pad=0.001) #,
    fig.savefig("%slife_inte_%s.png"%(figdir,suffix), bbox_inches='tight',pad_inches=0.08)

