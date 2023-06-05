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
behv   = ['remote','local']
numod= [chr(i) for i in range(97,115)]

def main_run():
    #draw_radiu_ts('%s/match_radiu_test.png'%figdir)
    draw_stacked_bar(6,'%s/match_remote_local_percent.pdf'%figdir)

def draw_stacked_bar(radiu,figname):
    titls = ['DJF','MAM','JJA','SON']
    nrow = 3 #6 #
    ncol = 2 #2 #
    bmlo = 0.3 #0.25 #
    
    ds = xr.open_dataset('%s/behv_season_6cyclone_6rad.nc'%outdir)
    total = ds['numb'].data[1,:,:,:].sum(axis=1) #[suffix,lev,behv,season]
    print(total[0,0])
    
    var = np.zeros([2,len(lev),len(lev),4], dtype=int)
    for nb in range(len(behv)):
        outfile = '%s/match_local_%s_%dcyclone_%drad.nc'%(
                outdir,behv[nb],radiu1,radiu)
        ds = xr.open_dataset(outfile)
        var[nb,:,:,:] = ds['numb'].data #[local_lev,remote_lev,season]

    colors = ['c','y','violet']
    x = np.arange(len(titls))
    fig = plt.figure(figsize=(12,10),dpi=300)
    ax = fig.subplots(nrow, ncol)
    for nb in range(len(behv)):
        for nl1 in range(len(lev)):
            axe = ax[nl1][nb]
            axe.set_title('(%s) %dhPa local genesis with %s'%(numod[2*nl1+nb],
                lev[nl1],behv[nb]),fontsize=title_font,fontdict=font)
            axe.set_ylabel('',fontsize=label_font,fontdict=font)
            for nl2 in range(len(lev)):
                axe.bar(x+0.3*nl2, var[nb,nl1,nl2,:], width=0.3, align='edge', 
                        bottom=0, color=colors[nl2])
                for nm,tota in enumerate(total[nl1,:]):
                    if round(var[nb,nl1,nl2,nm])>0:
                        axe.text(x[nm]+0.3*nl2+0.15, var[nb,nl1,nl2,nm]*2./3., '%d\n%d%%'%(
                            round(var[nb,nl1,nl2,nm]),round(var[nb,nl1,nl2,nm]*100/tota)),
                            ha='center', color='k', weight='bold', size=10)
            axe.set_xticks(x+0.45)
            axe.set_xticklabels(titls)
            axe.set_yticks([])
        
            if nb == 0 and nl1 == 0:
                patchs = []
                for nl2 in range(len(lev)):
                    patchs.append(Patch(color=colors[nl2],label=lev[nl2]))
                axe.legend(handles=patchs,handleheight=1.5,ncol=3)

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

def draw_radiu_ts(figname):
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
            axe.grid(True, which="both", axis='y',color='grey', linestyle='--', linewidth=1)
    axe.legend(titls)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()


