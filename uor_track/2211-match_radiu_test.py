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
path = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev'
outdir = '/home/ys17-23/Extension2/renql/project/uor_track/mdata'
figdir = '/home/ys17-23/Extension2/renql/project/uor_track/fig'
radiu1 = 6
radiu2 = [3,6,9,12]
suffix = ["%dlocal"%radiu1,"%doutside"%radiu1]
behv   = ['lysis','moveout']

def main_run():
    outfile = '%s/match_local_local_%dcyclone_%drad.nc'%(
            outdir,radiu1,radiu2)

def draw_stacked_bar(outfile,figname):
    titls = ['DJF','MAM','JJA','SON']
    nrow = 3 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    fig = plt.figure(figsize=(8,12),dpi=300)
    ax = fig.subplots(nrow, ncol)

    ds = xr.open_dataset('%s/behv_season_6cyclone_6rad.nc'%outdir)
    total = ds['numb'].data[0,:,:,:].sum(axis=1) #[suffix,lev,behv,season]
    print(total[0,0])
    
    ds = xr.open_dataset(outfile)
    var = ds['numb'].data #[local_lev,remote_lev,season]

    colors = ['c','y','violet']
    x = np.arange(len(titls))
    for nc in range(len(lev)):
        axe = ax[nc]
        axe.set_title('%dhPa local local %d'%(lev[nc],radiu2),fontsize=title_font,fontdict=font)
        axe.set_ylabel('number',fontsize=label_font,fontdict=font)
        
        for nl in range(len(lev)):
            axe.bar(x+0.3*nl, var[nc,nl,:], width=0.3, align='edge', 
                    bottom=0, color=colors[nl])
            
            for nm,tota in enumerate(total[nc,:]):
                if round(var[nc,nl,nm])>0:
                    axe.text(x[nm]+0.3*nl+0.15, var[nc,nl,nm]*2./3., '%d\n%d%%'%(
                        round(var[nc,nl,nm]),round(var[nc,nl,nm]*100/tota)),
                        ha='center', color='k', weight='bold', size=10)

        axe.set_xticks(x+0.45)
        axe.set_xticklabels(titls)
        
        if nc == 0:
            patchs = []
            for nb in range(len(lev)):
                patchs.append(Patch(color=colors[nb],label=lev[nb]))
            axe.legend(handles=patchs,handleheight=1.5,ncol=3)

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()


