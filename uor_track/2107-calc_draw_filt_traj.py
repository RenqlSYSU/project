#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
import sys, os, subprocess, linecache, gc
from datetime import datetime
from scipy import stats
from renql import monthly_calc, life_intensity, cyc_filter 
import matplotlib.pyplot as plt
from matplotlib import colors

if len(sys.argv) < 2 :
    flats = 45 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 85 #int(sys.argv[4])
    flonr = 110 #int(sys.argv[5])
    time = 48 # threshold, hour
    option = 2
    prefix  = "ff"
    suffix0 = "_5_2545-60110"
else:
    flats = int(sys.argv[1])
    flatn = int(sys.argv[2])
    flonl = int(sys.argv[3])
    flonr = int(sys.argv[4])
    time  = int(sys.argv[5])
    prefix = sys.argv[5]

suffix=str(option)+"_"+str(flatn)+str(flatn)+"-"+str(flonl)+str(flonr)
os.chdir("/home/users/qd201969/uor_track/fig")
path = '/home/users/qd201969/ERA5-1HR-lev/'
lev  = [850,500,250]
nday=[365,31,28,31,30,31,30,31,31,30,31,30,31]
month=["ALL","Jan","Feb","Mar","Apr","May","Jun",\
        "Jul","Aug","Sep","Oct","Nov","Dec"]
frae=744

numb = np.empty( [len(lev)],dtype=int )  
for nl in range(len(lev)):
    filname = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0
    
    if not os.path.isfile(filname+"_"+suffix) :
        numb[nl] = cyc_filter.line_filt(filname,flats,flatn,flonl,flonr,time=time)
    
    monthly_calc.draw_month_traj(frae, month[1:len(month)],nday[1:len(month)],
            filname+"_"+suffix,flats,flatn,flonl,flonr)
'''
com = "python /home/users/qd201969/uor_track/2107-draw_panel_traj.py %s %s_%s %d %d %d %d"\
        %(prefix,suffix0,suffix,flats,flatn,flonl,flonr)
ret = subprocess.Popen(com,shell=True)
ret.wait()
subprocess.run("rm %s%s_*_1980-2020%s_%s*"%(path,prefix,suffix0,suffix),shell=True)
'''

