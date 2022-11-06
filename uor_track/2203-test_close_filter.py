#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-6060
then judge whether the cyclone is NTN, STN or LYS in this region
if not, it is considered to be able to pass over the TP
plot the filter box in the trajectory figure

20211014
renql
'''

import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import sys, os, subprocess, linecache, gc
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import cartopy.crs as ccrs
from renql import cyc_filter, monthly_calc, life_intensity, add_field_filter

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

lev = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
path1 = '/home/users/qd201969/TRACK-1.5.2/'
figdir = "/home/users/qd201969/uor_track/fig/"
os.chdir(figdir)

prefix = "fftadd"
suffix = "_5_2545-60110_2_2545-6060"

def main_run():
    for nl in lev: 
        filname = "%s%s_%d_1980-2020%s"%(path,prefix,nl,suffix)
        numb = add_field_filter.filt_close_cyclone(filname,100)
        
        if numb > 0:
            draw_one_traj("%s_%dclose"%(filname,100))

def draw_one_traj(filname):
    ret=subprocess.Popen("%sutils/bin/tr2nc %s s %sutils/TR2NC/tr2nc.meta"%(
        path1,filname,path1),shell=True)
    ret.wait()
    f=cf.read(filname+'.nc')
    print(f)
    g = f[2]
    g = g*1e5
    print(g)
    #subprocess.run("rm %s.nc"%filname,shell=True) 
    
    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    term = filname.split("/")[-1]
    lev = term.split("_")[1]
    cfp.setvars(file="traj-%s.png"%term)
    cfp.gopen()
    cfp.mapset(lonmin=0, lonmax=150, latmin=15, latmax=70)
    
    cfp.levs(manual=[1500,3000,4500])
    cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
    cfp.levs()
    
    cfp.cscale('precip2_17lev')
    cfp.levs(min=0.0, max=8.0, step=0.5)
    cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
    cfp.traj(g, zorder=0, legend_lines=True, linewidth=2, 
            title='western close %s'%lev ,colorbar_title=' Relative Vorticity (Hz) * 1e5') 
    cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./traj-"+prefix+"_"+suffix+".png",shell=True) 

if __name__=='__main__':
    main_run()
