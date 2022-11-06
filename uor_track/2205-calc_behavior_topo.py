#!/usr/bin/env python
'''
1. read /home/users/qd201969/gtopo30_0.9x1.25.nc
    set the phis larger than 1500m as 1
2. if the 6degree cycle of cyclone center has 
    encounter 1, then set this cyclone as the one
    passing through the TP
3. Using location of the first point to define 
    local and remote cyclones

20220518
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
import matplotlib as mpl
import seaborn as sns

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev'
outdir = "/home/users/qd201969/uor_track/mdata"
figdir = "/home/users/qd201969/uor_track/fig"
prefix = 'ff'
radiu = 6
behv = ["%dtotal"%radiu,"%dlocal"%radiu,"%doutside"%radiu]

def main_run():
    thre = 1500
    varname = ['lifetime','distance','max vor','mean vor']
    for nint in [0,2]:
        draw_seasonal_box_1x2(varname[nint],nint)
    '''
    if not os.path.isfile('%s/tp_loca_%d.txt'%(outdir,thre)):
        write_tp_grid(thre,'%s/tp_loca_%d.txt'%(outdir,thre))
    
    for nl in lev:
        tp_filt_cyclone('%s/%s_%d_1980-2020'%(path,prefix,nl),
            '%s/tp_loca_%d.txt'%(outdir,thre))
    for nb in range(0,len(behv),1):
      # calc statistics and draw figure 
        com = "sh ~/uor_track/control_era5_1hr_track.sh %s 2 1 _%s"\
            %(prefix,behv[nb])
        ret=subprocess.Popen(com,shell=True)
        ret.wait()
    '''

def draw_seasonal_box_1x2(varname,nint):
    # nint: 0 lifetime, 1 distance, 2 max-vor
    # 3 mean-vor, 4 min-pres, 5 mean-pres
    nrow = 1 #6 #
    ncol = 2 #2 #
    bmlo = 0.35 #0.25 #
    titls = ['DJF','MAM','JJA','SON']
    
    fig = plt.figure(figsize=(12,5),dpi=150)
    ax = fig.subplots(nrow, ncol)
    
    for nc in range(2):
        p90 = np.zeros([len(lev),len(titls)], dtype=float ) # 4 or 12 
        dicts = {'lev':[],'season':[],varname:[]}
        for nl in range(len(lev)):
            filname  = '%s/fftadd_%d_1980-2020_%s'%(path,lev[nl],behv[nc+1])
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
        axe.set_title(behv[nc+1],fontsize=title_font,fontdict=font)
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

def get_x_coordinates_of_seaborn_boxplot(ax):
    display_coordinates = []
    inv = ax.transData.inverted()    
    for c in ax.get_children():
        if type(c) == mpl.patches.PathPatch:
            display_coordinates.append((c.get_extents().x0+c.get_extents().x1)/2)
    return inv.transform(tuple(display_coordinates)) 

def sea2mon(ns):
    if ns == 'DJF':
        lst = [12,1,2]
    if ns == 'MAM':
        lst = [3,4,5]
    if ns == 'JJA':
        lst = [6,7,8]
    if ns == 'SON':
        lst = [9,10,11]
    return lst

def tp_filt_cyclone(filname,tp_file):
    if os.path.exists(filname+"_"+behv[0]):
        print('%s_%s exists'%(filname,behv[0]))
        return
    else:
        print('handle %s_%s'%(filname,behv[0]))
    loca = np.loadtxt(tp_file,usecols = (0,1))
    #prefix = filname.split("_",1)[0].rsplit("/")[-1]

    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term0 = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term0[0]))
    number = [int(term0[0]),]
   
    # create output file and tid lists for different cyclones
    outfile=[]
    for nb in range(0,len(behv),1):
        outfile.append(open(filname+"_"+behv[nb],"w"))
        outfile[nb].write(line1)
        outfile[nb].write(line2)
        outfile[nb].write(line3)
        outfile[nb].write(line4)
        locals()['tid_%s'%behv[nb]]=[] 

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
                if prefix in ['ffadd','fftadd']:
                    data.append(list(map(float, line.strip().replace(" &","").split(" "))))
                else:
                    data.append(list(map(float, line.strip().split(" "))))
           
            signal=-10
            for nl in range(0,len(data)-1,1):
                dist = np.square(data[nl][2]-loca[:,0])+np.square(data[nl][1]-loca[:,1]) 
                if dist.min()<radiu*radiu:
                    if nl == 0 : # and dist.min()<=0.6
                        signal = 1 # local
                    else:
                        signal = 2 # outside
                    break
            if signal > 0 : 
                locals().get('tid_%s'%behv[0]).append(term[2])
                outfile[0].write(lineid)
                outfile[0].write(linenum)
                for nll in value:
                    outfile[0].write(nll)
                locals().get('tid_%s'%behv[signal]).append(term[2])
                outfile[signal].write(lineid)
                outfile[signal].write(linenum)
                for nll in value:
                    outfile[signal].write(nll)

        line = ff.readline()

    ff.close()
    for nb in range(0,len(behv),1):
        number.append(len(locals().get('tid_%s'%behv[nb])))
        outfile[nb].seek(0,0) # Go back to the beginning of the file
        outfile[nb].write(line1)
        outfile[nb].write(line2)
        outfile[nb].write(line3)
        outfile[nb].write("TRACK_NUM %9d %s"%(number[nb+1],term0[1]))

        print("%s : %d" %(outfile[nb].name,number[nb+1]))
        outfile[nb].close()
    
    return number

def write_tp_grid(thre,figdir):
    lats = 23 #int(sys.argv[2])
    latn = 45 #int(sys.argv[3])
    lonl = 60 #int(sys.argv[4])
    lonr = 110 #int(sys.argv[5])
    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat).load()
    phis = phis/9.8 # transfer from m2/s2 to m
    
    dim = phis.shape
    f = open(figdir,'w')
    for i in range(0,dim[0]):
        for j in range(0,dim[1]):
            if phis[i,j]>thre:
                f.write("%.2f %.2f %.2f \n"%(
                    phis[i,j].lat,phis[i,j].lon,phis[i,j]))
    f.close()

if __name__=='__main__':
    main_run()

