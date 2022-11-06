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
import seaborn as sns
from datetime import datetime, timedelta
from renql import life_intensity
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
outdir = '/home/ys17-23/Extension2/renql/uor_track/mdata'
figdir = '/home/ys17-23/Extension2/renql/uor_track/fig'
radiu1 = 6
radiu2 = 6 # distance match 
suffix = ["%dlocal"%radiu1,"%doutside"%radiu1]

def main_run():
    varname = ['lifetime','distance','max vor','mean vor']
    for nint in [0,2]:
        draw_seasonal_box_3x1(varname[nint],nint)
    ''' 
    for nl1 in lev:
        for nl2 in lev:
            outfile = '%s/fftadd_match_%dlocal_%dremote_%ddist'%(
                outdir,nl1,nl2,radiu2)
            match_local_remote('%s/fftadd_%d_1980-2020_%dlocal'%(path,nl1,radiu1),
                '%s/fftadd_%d_1980-2020_%doutside'%(path,nl2,radiu1),radiu2,outfile)
    
    outfile = '%s/match_local_season_%dcyclone_%drad.nc'%(
            outdir,radiu1,radiu2)
    calc_month(outfile)
    draw_stacked_bar(outfile,'%s/bar_match_local_%dcyc_%drad.png'%(
        figdir,radiu1,radiu2))
    '''
def match_local_remote(filname1,filname2,dist,outfilename):
    if os.path.exists(outfilename):
        print('%s exists'%outfilename)
        return
    else:
        print('handle %s'%outfilename)
    
    ff = open(filname1,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    outfile = open(outfilename,"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write(line4)

    log = open('%s_log'%outfilename,'w')
    var = read_cyclone(filname2)
    tid=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            value=[]
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
            data = list(map(float, value[0].strip().replace(" &","").split(" ")))
           
            signal = match_cyclone2(data[0],data[1],data[2],var,dist) 
            if len(signal) >= 1:
                tid.append(term[2])
                outfile.write(lineid)
                outfile.write(linenum)
                for nll in value:
                    outfile.write(nll)
                log.write('%d %d %d %d %s\n'%(data[0],data[1],data[2],len(signal),signal))

        line = ff.readline()
    ff.close()
    log.close()
    if len(tid) == 0:
        outfile.close()
        subprocess.run('rm %s*'%outfilename,shell=True)
    else:
        outfile.seek(0,0) # Go back to the beginning of the file
        outfile.write(line1)
        outfile.write(line2)
        outfile.write(line3)
        outfile.write("TRACK_NUM %9d ADD_FLD    0   0 &" %len(tid))
        print("filter cyclone number in %s : %d" %(outfile.name,len(tid)))
        outfile.close()
    return len(tid)

def read_cyclone(filname):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    tid=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float,line.strip().replace(" &","").split(" ")))
                tid.append([data[0],data[1],data[2],nl+1,num])
        line = ff.readline()
    ff.close()
    tid = np.array(tid)
    print('%s %s'%(filname,str(tid.shape)))
    return tid

def match_cyclone2(time,lon,lat,var,dist):
    tid = []
    ind = np.argwhere(var[:,0]==time)[:,0]
    if len(ind)>0:
        print(len(ind))
        term = var[ind,:]
        distance = np.square(term[:,2]-lat)+np.square(term[:,1]-lon)
        print(distance)
        ind2 = np.argwhere(distance<=dist*dist)[:,0]
        if len(ind2)>0:
            tid = ['%d/%d'%(term[i,3],term[i,4]) for i in ind2]
    print('%f %f %f : match number %d %s'%(time,lon,lat,len(tid),tid))
    return tid 

def calc_month(outfile):
    if os.path.exists(outfile):
        print('%s exists'%outfile)
        return
    else:
        print('handle %s'%outfile)

    var = np.zeros([len(lev),len(lev),4], dtype=int)
    for nl1 in range(len(lev)):
        for nl2 in range(len(lev)):
            filname = '%s/fftadd_match_%dlocal_%dremote_%ddist'%(
                outdir,lev[nl1],lev[nl2],radiu2)
            
            if os.path.exists(filname):
                var[nl1,nl2,:] = calc_season_cyclone(filname) 
            else:
                var[nl1,nl2,:] = np.zeros([4])

    var = xr.DataArray(var/41.0)
    print(var)
    ds1 = var.to_dataset(name="numb")
    ds1.to_netcdf(outfile,'w')

def calc_season_cyclone(filname):
    var = np.zeros([4], dtype=int ) # 4 or 12 

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term0 = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term0[0]))
    
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
            if data[-1][1] > data[0][1]:
                if sum(i.month==ct1[0].month for i in ct1)/len(ct1) >= 0.5 : 
                    var[mon2sea(ct1[0].month)] += 1
                else:
                    var[mon2sea(ct1[-1].month)] += 1
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
        axe.set_title('%dhPa local'%lev[nc],fontsize=title_font,fontdict=font)
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

def draw_seasonal_box_3x1(varname,nint):
    # nint: 0 lifetime, 1 distance, 2 max-vor
    # 3 mean-vor, 4 min-pres, 5 mean-pres
    nrow = 3 #6 #
    ncol = 1 #2 #
    bmlo = 0.35 #0.25 #
    titls = ['DJF','MAM','JJA','SON']
    
    fig = plt.figure(figsize=(6,12),dpi=150)
    ax = fig.subplots(nrow, ncol)
    
    for nc in range(len(lev)):
        p90 = np.zeros([len(lev),len(titls)], dtype=float ) # 4 or 12 
        dicts = {'lev':[],'season':[],varname:[]}
        for nl in range(len(lev)):
            filname = '%s/fftadd_match_%dlocal_%dremote_%ddist'%(
                outdir,lev[nc],lev[nl],radiu2)
            var = life_intensity.calc_one_variable(filname,nint,
                flats=10,flatn=60,flonl=50,flonr=120)
            #var = life_intensity.calc_one_variable(filname,nint,
            #    flats=0,flatn=90,flonl=0,flonr=180)
            for nm in range(len(titls)):
                p90[nl,nm] = np.percentile(np.array(var[nm]),90)
                dicts['lev']  = dicts['lev']+[lev[nl]]*len(var[nm])
                dicts['season'] = dicts['season']+[titls[nm]]*len(var[nm])
                dicts[varname] = dicts[varname]+var[nm]
                print('%s %d %d 90th %f'%(titls[nm],lev[nl],len(var[nm]),p90[nl,nm]))
        df = pd.DataFrame(dicts)
        print(df)

        axe = ax[nc]
        axe.set_title('%dhPa local'%lev[nc],fontsize=title_font,fontdict=font)
        sns.boxplot(x='season', y=varname, hue='lev',hue_order=lev, 
            data=df, palette="Set1",ax=axe, showfliers=False, whis=0,
            showmeans=True,meanprops={"markerfacecolor":"k", "markeredgecolor":"k"})
        #axe.set_ylim(1,10)
        axe.set_ylabel(varname,fontsize=label_font,fontdict=font)
        axe.set_xlabel('',fontsize=label_font,fontdict=font)
        
        width = 0.8
        xloc = np.repeat(np.atleast_2d(np.arange(4)),3,axis=0
            )+np.array([[-1*width/3.0],[0],[width/3.0]])
        col = ['ro','bo','go']
        for nl in range(len(lev)):
            axe.plot(xloc[nl,:],p90[nl,:],col[nl])
            #axe.plot(xloc.flatten(),p90.flatten(),'ko')
        
    plt.legend([],[], frameon=False)
    plt.tight_layout()
    plt.savefig('%s/box_%s.png'%(figdir,varname), 
        bbox_inches='tight',pad_inches=0.01)

if __name__=='__main__':
    main_run()

