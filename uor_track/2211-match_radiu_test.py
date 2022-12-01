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
radiu2 = range(2,13,1) #[3,6,9,12]
suffix = ["%dlocal"%radiu1,"%doutside"%radiu1]
behv   = ['local','remote']

def main_run():
    draw_stacked_bar('%s/match_radiu_test.png'%figdir)

def draw_stacked_bar(figname):
    titls = ['DJF','MAM','JJA','SON']
    nrow = 3 #6 #
    ncol = 2 #2 #
    bmlo = 0.35 #0.25 #
    
    ds = xr.open_dataset('%s/behv_season_6cyclone_6rad.nc'%outdir)
    total = ds['numb'].data[0,:,:,:].sum(axis=1) #[suffix,lev,behv,season]
    print(total[0,0]) #[lev,season]
    
    var = np.zeros([2,len(radiu2),len(lev),4], dtype=int)
    for nb in range(len(behv)):
        for nr in range(len(radiu2)):
            outfile = '%s/match_local_%s_%dcyclone_%drad.nc'%(
                    outdir,behv[nb],radiu1,radiu2[nr])
            ds = xr.open_dataset(outfile)
            numb = ds['numb'].data #[local_lev,remote_lev,season]
            for nl,nll in enumerate([1,2,1]):
                var[nb,nr,nl,:] = numb[nl,nll,:]*100.0/total[nl,:]

    fig = plt.figure(figsize=(8,12),dpi=150)
    ax = fig.subplots(nrow, ncol)
    for nb in range(len(behv)):
        for nl in range(len(lev)):
            axe = ax[nl][nb]
            axe.set_title('%dhPa %s'%(lev[nl],behv[nb]),fontsize=title_font,fontdict=font)
            axe.set_ylabel('percent',fontsize=label_font,fontdict=font)
            for ns in range(4):
                #axe.plot(radiu2,var[nb,:,nl,ns],linewidth=3,marker="o",markersize=6)
                axe.plot([x**2 for x in radiu2],var[nb,:,nl,ns],linewidth=3,marker="o",markersize=6)
    axe.legend(titls)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()


