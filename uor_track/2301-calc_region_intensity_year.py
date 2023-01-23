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
radiu1 = 6
radiu2 = 6
nmonth = 12
behv   = ['upstm','TP','downstm']
lon = [40,65,105,130]
#lat = [15,60]
lat = [20,60]
nyear = 41
numod= [chr(i) for i in range(97,115)]

def main_run():
    outfile = '%s/region_intensity_month.nc'%(outdir)
    calc_region_intensity(outfile)
    draw_ts_3x1(outfile,'%s/intensity_%dcyc_%drad.png'%(
        figdir,radiu1,radiu2))

def calc_region_intensity(outfile):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    var = np.zeros([len(lev),len(behv),nyear,nmonth,2], dtype=float)
    for nl in range(len(lev)):
        filname = '%s/fftadd_%d_1980-2020_%dtotal'%(path,lev[nl],radiu1)
        var[nl,:,:,:,:] = calc_season_cyclone(filname,var[nl,:,:,:,:]) 
   
    var[:,:,:,:,0] = var[:,:,:,:,0]/np.where(
        var[:,:,:,:,1]==0, 1, var[:,:,:,:,1])
    var = xr.DataArray(var)
    print(var)
    ds1 = var.to_dataset(name="inte")
    ds1.to_netcdf(outfile,'w')

def calc_season_cyclone(filname,var):
    loca = np.loadtxt(tp_file,usecols = (0,1))

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
            
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float,line.strip().replace(" &","").split(" ")))
                ct1 = datetime.strptime(str(int(data[0])),'%Y%m%d%H')
                for signal in range(len(behv)):
                    if data[1]>lon[signal] and data[1]<lon[signal+1] and data[2]>lat[0] and data[2]<lat[1] :
                        var[signal,ct1.year-1981,ct1.month-1,0] += data[3]
                        var[signal,ct1.year-1981,ct1.month-1,1] += 1

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

def draw_ts_3x1(outfile,figname):
    nrow = 3 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    ds = xr.open_dataset(outfile)
    var = ds['inte'].data[:,:,:,:,0] 
    
    pcolor  = ["r","b",'g'] # change with option
    x = np.arange(1,nmonth+1,1)
    titls = ["Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec"]
        
    fig = plt.figure(figsize=(8,12),dpi=300)
    ax = fig.subplots(nrow, ncol)

    for nl in range(len(lev)):
        axe = ax[nl]
        axe.set_title('(%s) %dhPa'%(numod[nl],lev[nl]),
            fontsize=title_font,fontdict=font)
        axe.set_ylabel('intensity',fontsize=label_font,fontdict=font)
        for nb in range(len(behv)):
            axe.plot(x,var[nl,nb,:,:].mean(0),linewidth=2,marker="o",markersize=6,
                color=pcolor[nb])
            #axe.plot(x,var[nl,nb,:,:].mean(0)+np.std(var[nl,nb,:,:],0),linewidth=0.5,color=pcolor[nb])
            #axe.plot(x,var[nl,nb,:,:].mean(0)-np.std(var[nl,nb,:,:],0),linewidth=0.5,color=pcolor[nb])
        axe.grid(True, which="both", axis='y',color='grey', linestyle='--', linewidth=1)
        axe.set_xticks(x)
        axe.set_xticklabels(titls)
        
        ax2 = axe.twinx()
        tvalue,pvalue = stats.ttest_ind(var[nl,0,:,:],var[nl,1,:,:], axis=0, equal_var=False)
        ax2.plot(x,pvalue,linewidth=2,color='brown')
        tvalue,pvalue = stats.ttest_ind(var[nl,2,:,:],var[nl,1,:,:], axis=0, equal_var=False)
        ax2.plot(x,pvalue,linewidth=2,color='c')
        #ax2.plot([1,12],[0.01,0.01],linewidth=2,color='k')

    axe.legend(behv,loc='upper right')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname,bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

