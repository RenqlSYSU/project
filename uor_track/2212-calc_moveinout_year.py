#!/usr/bin/env python
import sys
import subprocess, os
import xarray as xr
import numpy as np
import pandas as pd 
from scipy import stats
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from datetime import datetime, timedelta
#matplotlib.use('Agg')

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

lev = [850, 500, 250]
path  = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev'
path1 = '/home/ys17-23/Extension2/renql/project/uor_track'
outdir = '/home/ys17-23/Extension2/renql/project/uor_track/mdata'
figdir = '/home/ys17-23/Extension2/renql/project/uor_track/fig'
tp_file='%s/tp_loca_1500.txt'%(outdir)
nmonth = 12
radiu1 = 3
radiu2 = 6
behv   = ['in','out']
slon = 60
elon = 110
nyear = 41
numod= [chr(i) for i in range(97,115)]

def main_run():
    #outfile = '%s/movein_moveout_month_%dcyclone_%drad_%d.nc'%(
    #        outdir,radiu1,radiu2,elon)
    outfile = '%s/movein_moveout_month_%dcyclone_%drad_year.nc'%(
            outdir,radiu1,radiu2)
    calc_lysis_percent(outfile)
    draw_ts_3x1(outfile,'%s/moveinout_%dcyc_%drad.png'%(
        figdir,radiu1,radiu2))
    #draw_6ts(outfile,'%s/moveinout_%dcyc_%drad.png'%(
    #    figdir,radiu1,radiu2))

def calc_lysis_percent(outfile):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    var = np.zeros([len(lev),len(behv),nyear,nmonth], dtype=int)
    for nl in range(len(lev)):
        filname = '%s/fftadd_%d_1980-2020_%dtotal'%(path,lev[nl],radiu1)
        var[nl,:,:,:] = calc_season_cyclone(filname) 

    var = xr.DataArray(var)
    print(var)
    ds1 = var.to_dataset(name="numb")
    ds1.to_netcdf(outfile,'w')

def calc_season_cyclone(filname):
    loca = np.loadtxt(tp_file,usecols = (0,1))
    var = np.zeros([len(behv),nyear,nmonth], dtype=int ) # 4 or 12 

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term0 = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term0[0]))
   
    prefix = filname.split("_",1)[0].rsplit("/")[-1]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            ct1=[]
            data=[]
            for nl in range(0,num,1):
                line = ff.readline()
                data.append(list(map(float,line.strip().replace(" &","").split(" "))))
                ct1.append(datetime.strptime(str(int(data[-1][0])),'%Y%m%d%H'))

            signal=-10
            if data[-1][1] > data[0][1] or data[0][1]>180 :
                dist = np.square(data[0][2]-loca[:,0])+np.square(data[0][1]-loca[:,1]) 
                if dist.min()>radiu2*radiu2:
                    signal = 0 # movein 
                    if sum(i.month==ct1[0].month for i in ct1)/len(ct1) >= 0.5 : 
                        var[signal,ct1[0].year-1980,ct1[0].month-1] += 1
                    else:
                        var[signal,ct1[-1].year-1980,ct1[-1].month-1] += 1

                dist = np.square(data[-1][2]-loca[:,0])+np.square(data[-1][1]-loca[:,1]) 
                if dist.min()>radiu2*radiu2 :#and data[-1][1]>115: 
                    signal = 1 # moveout
                    if sum(i.month==ct1[0].month for i in ct1)/len(ct1) >= 0.5 : 
                        var[signal,ct1[0].year-1980,ct1[0].month-1] += 1
                    else:
                        var[signal,ct1[-1].year-1980,ct1[-1].month-1] += 1

        line = ff.readline()
    ff.close()
    return var 

def mon2sea(nm):
# convert month to season index
    if nm in [1,2,12]:
        sea = 0
    elif nm in [3,4,5]:
        sea = 1
    elif nm in [6,7,8]:
        sea = 2
    elif nm in [9,10,11]:
        sea = 3
    return sea

def draw_6ts(outfile,figname):
    nrow = 1 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    ds = xr.open_dataset(outfile)
    var = ds['numb'].data #[lev,behv,season]
    
    plinsty = ["solid",(0,(4,2)),(0,(1,2))] # change with lev
    pcolor  = ["k","r","b"] # change with option
    x = np.arange(1,nmonth+1,1)
        
    fig = plt.figure(figsize=(12,8),dpi=300)
    axe = fig.subplots(nrow, ncol)
    axe.set_title('',fontsize=title_font,fontdict=font)
    axe.set_ylabel('number',fontsize=label_font,fontdict=font)
    axe.set_xlabel('month',fontsize=label_font,fontdict=font)
    for nb in range(len(behv)):
        for nl in range(len(lev)):
            axe.plot(x,var[nl,nb,:],linewidth=3,marker="o",markersize=6,
                color=pcolor[nl],linestyle=plinsty[nb],
                label='%d%s'%(lev[nl],behv[nb]))
    axe.grid(True, which="both", axis='x',color='grey', linestyle='--', linewidth=1)
    axe.legend(ncol=2,fontsize=label_font)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

def draw_ts_3x1(outfile,figname):
    nrow = 3 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    ds = xr.open_dataset(outfile)
    var = ds['numb'].data #[lev,behv,season]
    
    pcolor  = ["r","b"] # change with option
    x = np.arange(1,nmonth+1,1)
    titls = ["Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec"]
        
    fig = plt.figure(figsize=(8,12),dpi=300)
    ax = fig.subplots(nrow, ncol)

    for nl in range(len(lev)):
        axe = ax[nl]
        axe.set_title('(%s) %dhPa'%(numod[nl],lev[nl]),
            fontsize=title_font,fontdict=font)
        axe.set_ylabel('number',fontsize=label_font,fontdict=font)
        for nb in range(len(behv)):
            axe.plot(x,var[nl,nb,:,:].mean(0),linewidth=2,marker="o",markersize=6,
                color=pcolor[nb])
            #axe.plot(x,var[nl,nb,:,:].mean(0)+np.std(var[nl,nb,:,:],0),linewidth=0.5,color=pcolor[nb])
            #axe.plot(x,var[nl,nb,:,:].mean(0)-np.std(var[nl,nb,:,:],0),linewidth=0.5,color=pcolor[nb])
        axe.grid(True, which="both", axis='y',color='grey', linestyle='--', linewidth=1)
        axe.set_xticks(x)
        axe.set_xticklabels(titls)
        
        tvalue,pvalue = stats.ttest_ind(var[nl,0,:,:],var[nl,1,:,:], axis=0, equal_var=False)
        ax2 = axe.twinx()
        ax2.set_ylabel('p-value',fontsize=label_font,fontdict=font)
        ax2.plot(x,pvalue,linewidth=2,color='gray')
        #ax2.plot([1,12],[0.01,0.01],linewidth=2,color='k')

    axe.legend(behv,loc='lower right')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

