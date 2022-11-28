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
path  = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev'
path1 = '/home/ys17-23/Extension2/renql/project/uor_track'
outdir = '/home/ys17-23/Extension2/renql/project/uor_track/mdata'
figdir = '/home/ys17-23/Extension2/renql/project/uor_track/fig'
tp_file='%s/tp_loca_1500.txt'%(outdir)
radiu1 = 6
radiu2 = 6 # use to filte lysis cyclone
suffix = ["%dlocal"%radiu1,"%doutside"%radiu1]
behv   = ['lysis','moveout']

def main_run():
    '''
    outfile = '%s/behv_season_%dcyclone_%drad.nc'%(
            outdir,radiu1,radiu2)
    calc_lysis_percent(outfile)
    draw_stacked_bar(outfile,'%s/bar_lysis_percent_%dcyc_%drad.png'%(
        figdir,radiu1,radiu2))
    '''
    for nc in range(len(suffix)):
        for nl in range(len(lev)):
            filname = '%s/ff_%d_1980-2020_%s'%(path,lev[nl],suffix[nc])
            write_moveout_cyclone(filname,1)
      # calc statistics and draw figure 
        com = "bash %s/control_era5_1hr_track.sh ff 2 1 _%s_moveout"\
            %(path1,suffix[nc])
        ret=subprocess.Popen(com,shell=True)
        ret.wait()

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
            var[nc,nl,:,:] = calc_season_cyclone(filname) 

    var = xr.DataArray(var/41.0)
    print(var)
    ds1 = var.to_dataset(name="numb")
    ds1.to_netcdf(outfile,'w')

def calc_season_cyclone(filname):
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
            
            ct1=[]
            data=[]
            for nl in range(0,num,1):
                line = ff.readline()
                data.append(list(map(float,line.strip().split(" "))))
                ct1.append(start+timedelta(hours=int(data[-1][0])))

            signal=-10
            if data[-1][1] > data[0][1]:
                dist = np.square(data[-1][2]-loca[:,0])+np.square(data[-1][1]-loca[:,1]) 
                if dist.min()<radiu2*radiu2:
                    signal = 0 # lysis
                else:
                    signal = 1 # moveout
                    
                if sum(i.year==1996 for i in ct1)/len(ct1) >= 0.5 : 
                    if sum(i.month==ct1[0].month for i in ct1)/len(ct1) >= 0.5 : 
                        var[signal,mon2sea(ct1[0].month)] += 1
                    else:
                        var[signal,mon2sea(ct1[-1].month)] += 1
        line = ff.readline()
    ff.close()
    return var 

def write_moveout_cyclone(filname,nb):
    if os.path.exists(filname+"_"+behv[nb]):
        print('%s_%s exists'%(filname,behv[nb]))
        return
    else:
        print('handle %s_%s'%(filname,behv[nb]))
    loca = np.loadtxt(tp_file,usecols = (0,1))

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term0 = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term0[0]))
   
    outfile = open(filname+"_"+behv[nb],"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write(line4)

    tid=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data=[]
            value=[]
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
                data.append(list(map(float, line.strip().split(" "))))
           
            signal=-10
            if data[-1][1] > data[0][1]:
                dist = np.square(data[-1][2]-loca[:,0])+np.square(data[-1][1]-loca[:,1]) 
                if dist.min()<radiu2*radiu2:
                    signal = 0 # lysis
                else:
                    signal = 1 # moveout
                    
            if signal == nb : 
                tid.append(term[2])
                outfile.write(lineid)
                outfile.write(linenum)
                for nll in value:
                    outfile.write(nll)

        line = ff.readline()

    ff.close()
    outfile.seek(0,0) # Go back to the beginning of the file
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write("TRACK_NUM %9d %s"%(len(tid),term0[1]))

    print("%s : %d" %(outfile.name,len(tid)))
    outfile.close()
    

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

def draw_stacked_bar(outfile,figname):
    titls = ['DJF','MAM','JJA','SON']
    nrow = 2 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    
    fig = plt.figure(figsize=(8,12),dpi=300)
    ax = fig.subplots(nrow, ncol)

    ds = xr.open_dataset(outfile)
    var = ds['numb'].data #[suffix,lev,behv,season]
    total = np.sum(var,axis=2)

    colors = ['c','y']
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
                    axe.text(x[nm]+0.3*nl+0.15, bot+var[nc,nl,nb,nm]+yoffset2[nc], 
                        '%d%%'%(round(var[nc,nl,nb,nm]*100/tota)),
                        #round(var[nc,nl,nb,nm]),
                        ha='center', color='k', weight='bold', size=10)
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
            axe.legend(handles=patchs,handleheight=1.5,ncol=2)

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figname, 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

