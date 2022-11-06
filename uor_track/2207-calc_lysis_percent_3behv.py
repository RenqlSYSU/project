#!/usr/bin/env python
import sys
import subprocess, os
import xarray as xr
import numpy as np
import pandas as pd 
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
path = '/home/users/qd201969/ERA5-1HR-lev'
outdir = '/home/users/qd201969/uor_track/mdata'
figdir = '/home/users/qd201969/uor_track/fig'
tp_file='%s/tp_loca_1500.txt'%(outdir)
radiu1 = 6
radiu2 = 6 # use to filte lysis cyclone
suffix = ["%dlocal"%radiu1,"%doutside"%radiu1]
behv   = ['lysis','EA','moveout']
lonxs  = [115,110]
time = 48

def main_run():
    outfile = '%s/behv3_%dcyclone_%drad_%d%d_%d.nc'%(
            outdir,radiu1,radiu2,lonxs[0],lonxs[1],time)
    calc_lysis_percent(outfile)
    draw_stacked_bar(outfile,'%s/bar_lysis_percent_%dcyc_%drad.png'%(
        figdir,radiu1,radiu2),True)

def calc_lysis_percent(outfile):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    var = np.zeros([len(suffix),len(lev),len(behv),4], dtype=int)
    for nc in range(len(suffix)):
        for nl in range(len(lev)):
            filname = '%s/ff_%d_1980-2020_%s'%(path,lev[nl],suffix[nc])
            var[nc,nl,:,:] = calc_season_cyclone(filname,lonxs[nc]) 

    var = xr.DataArray(var/41.0)
    print(var)
    ds1 = var.to_dataset(name="numb")
    ds1.to_netcdf(outfile,'w')

def calc_season_cyclone(filname,lonx):
    loca = np.loadtxt(tp_file,usecols = (0,1))
    var = np.zeros([len(behv),4], dtype=int ) # 4 or 12 

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term0 = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term0[0]))
   
    start=datetime(1995, 11, 30, 23, 00)
    prefix = filname.split("_",1)[0].rsplit("/")[-1]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
           
            if num >= time:
                ct1=[]
                for nl in range(0,num,1):
                    line = ff.readline()
                    data = list(map(float,line.strip().split(" ")))
                    ct1.append(start+timedelta(hours=int(data[0])))

                signal=-10
                dist = np.square(data[2]-loca[:,0])+np.square(data[1]-loca[:,1]) 
                if dist.min()<radiu2*radiu2:
                    signal = 0 # lysis
                elif data[1]>lonx:
                    signal = 1
                else:
                    signal = 2 # moveout
                    
                if sum(i.year==1996 for i in ct1)/len(ct1) >= 0.5 : 
                    if sum(i.month==ct1[0].month for i in ct1)/len(ct1) >= 0.5 : 
                        var[signal,mon2sea(ct1[0].month)] += 1
                    else:
                        var[signal,mon2sea(ct1[-1].month)] += 1
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

def draw_stacked_bar(outfile,figname,percent):
    titls = ['DJF','MAM','JJA','SON']
    nrow = 2 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    fig = plt.figure(figsize=(8,12),dpi=300)
    ax = fig.subplots(nrow, ncol)

    ds = xr.open_dataset(outfile)
    var = ds['numb'].data #[suffix,lev,behv,season]
    total = np.sum(var,axis=2)

    colors = ['c','y','orange']
    dash   = ['--','//','x']
    x = np.arange(len(titls))
    yoffset1 = [1,0.8]
    yoffset2 = [-10,-5]
    for nc in range(len(var)):
        axe = ax[nc]
        axe.set_title(suffix[nc],fontsize=title_font,fontdict=font)
        axe.set_ylabel('number',fontsize=label_font,fontdict=font)
        
        for nl in range(len(lev)):
            bottom = np.zeros([4])
            for nb in range(len(behv)):
                axe.bar(x+0.3*nl, var[nc,nl,nb,:], width=0.3, align='edge', 
                        bottom=bottom, color=colors[nb])#, hatch=dash[nl]
                bottom += var[nc,nl,nb,:]
            
            for nm,tota in enumerate(total[nc,nl,:]): 
                axe.text(x[nm]+0.3*nl+0.15, tota+yoffset1[nc], round(tota),
                      ha='center', color='k', weight='bold', size=10)

                bot = 0
                for nb in range(len(behv)):
                    if percent :
                        txt = '%d%%'%(round(var[nc,nl,nb,nm]*100/tota))
                    else:
                        txt = round(var[nc,nl,nb,nm]) 
                    if var[nc,nl,nb,nm] > abs(yoffset2[nc]):
                        axe.text(x[nm]+0.3*nl+0.15, bot+var[nc,nl,nb,nm]+yoffset2[nc], 
                            txt, ha='center', color='k', weight='bold', size=10)
                    bot += var[nc,nl,nb,nm]

        axe.set_xticks(x+0.45)
        axe.set_xticklabels(titls)
        
#        y_offset = -3
#        for bar in axe.patches:
#            axe.text(bar.get_x()+bar.get_width()/2,
#                  bar.get_height()+bar.get_y()+y_offset,round(bar.get_height()),
#                  ha='center', color='k', weight='bold', size=10)
#                  #backgroundcolor='w')

        if nc == 0:
            patchs = []
            #for nl in range(len(lev)):
            #    patchs.append(Patch(facecolor='w',hatch=dash[nl],label='%dhPa'%lev[nl]))
            for nb in range(len(behv)):
                patchs.append(Patch(color=colors[nb],label=behv[nb]))
            axe.legend(handles=patchs,handleheight=1.5,ncol=len(behv))

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

