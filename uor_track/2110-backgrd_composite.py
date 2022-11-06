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
import pandas as pd
import sys, os, subprocess, linecache, gc
from datetime import datetime

def composite_time(filname,flats,flatn,flonl,flonr,alltime):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    ct1=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    a = str(int(data[0]))
                    if int(a[-2:]) in range(0,6,1):
                        b='00'
                    elif int(a[-2:]) in range(6,12,1):
                        b='06'
                    elif int(a[-2:]) in range(12,18,1):
                        b='12'
                    elif int(a[-2:]) in range(18,24,1):
                        b='18'
                    a = a[:-2]+b
                    ct1.append(datetime.strptime(a,'%Y%m%d%H'))
        line = ff.readline()
    ff.close()
    ct2=list(set(ct1))
    # print total, unrepeated, repeated, repeated/total, unrepeated/alltime
    print("total, unrepeated, repeated, repeated/total, unrepeated/alltime")
    if len(ct1) == 0:
        print("all time when cyclone in %d-%dE,%d-%dN : %d"\
                %(flonl,flonr,flats,flatn,len(ct1)))
    else:
        print("all time when cyclone in %d-%dE,%d-%dN : %d, %d, %d, %.2f%%, %.2f%%"\
                %(flonl,flonr,flats,flatn,len(ct1),len(ct2),(len(ct1)-len(ct2)),\
                (len(ct1)-len(ct2))*100/len(ct1),len(ct2)*100/alltime))
    print("")
    return ct2

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 27 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 90 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "fft"
    season = 0 # 0 monthly, 1 seasonal
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
    season = int(sys.argv[7])
    time = int(sys.argv[8])

suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonl)
figdir = "/home/users/qd201969/uor_track/fig/behv3_month_%dh_%s"%(time,suffix)
fileout="/home/users/qd201969/uor_track/mdata/comp_6h_"

behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
levc = [850,500,250,200]
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/"#t/ERA5_NH_t_1989.nc

ftime  = pd.date_range(start='1979-12-01 00',end='2021-01-31 23', freq='6H',closed=None)
alltime= len(ftime)

for nl in range(0,len(lev),1):
    for nr in range(1,len(behv),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix+"_"+behv[nr]
        locals()['ct_%d_%s'%(lev[nl],behv[nr])]=composite_time(filname,flats,flatn,flonl,flonr,alltime)

lonl=0  #0  #
lonr=150#360#
lats=0  #
latn=90 #
ds  = xr.open_dataset(datapath+"t/ERA5_NH_t_1989.nc")
lat = ds.latitude
lon = ds.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
for varname in ['u','v','t','z']:
    var  = np.empty( [len(lev),len(behv),len(levc),len(ilat),len(ilon)],dtype=float )  
    numb = np.empty( [len(lev),len(behv)],dtype=int ) 
    for nl in range(0,len(lev),1):
        for nr in range(0,len(behv),1):
            term = np.zeros( [len(levc),len(ilat),len(ilon)],dtype=float ) 
            print("handle %d %s"%(lev[nl],behv[nr]))
            if nr == 0:
                for year in range(1979,2021,1):
                    ds  = xr.open_dataset("%s%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,year))
                    term= term + ds[varname].sel(level=levc,longitude=ilon,latitude=ilat).sum('time')
                allt = len(pd.date_range(start='1979-01-01 00',end='2020-12-31 23', freq='6H',closed=None))
                var[nl,nr,:,:,:] = term/allt
                numb[nl,nr]=allt
                del term
            else:
                ct = locals().get('ct_%d_%s'%(lev[nl],behv[nr]))
                if len(ct) == 0:
                    numb[nl,nr]=0
                    continue
                ctda = xr.DataArray(ct, coords=[ct], dims=["time"])
                for year in range(1979,2021,1):
                    ctda1 = ctda.sel(time=ctda.time.dt.year.isin(year))
                    #print("%d numb %d"%(year,len(ctda1)))
                    if len(ctda1) == 0:
                        continue
                    ds  = xr.open_dataset("%s%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,year))
                    term= term + ds[varname].sel(time=ctda1,level=levc,longitude=ilon,latitude=ilat).sum('time')
                var[nl,nr,:,:,:] = term/len(ctda.sel(time=ctda.time.dt.year.isin(range(1979,2021,1))))
                numb[nl,nr]=len(ctda.sel(time=ctda.time.dt.year.isin(range(1979,2021,1))))
                del ct, ctda, term
        
    ds = xr.Dataset(
            {
                "numb": (["lev3", "behv"], numb),
                "var" : (["lev3", "behv", "level", "lat","lon"], var),
                },
            coords={
                "lev3" : (["lev3"],lev),
                "behv" : (["behv"],behv),
                "level": (["level"],levc),
                "lat"  : (["lat"],ilat.data),
                "lon"  : (["lon"],ilon.data),
                },
            )
    ds.attrs["description"] = "composite background for diff behavior of 3levs cyclones"
    ds.to_netcdf(fileout+varname+".nc","w")
    del ds,var,numb
        
print("")
print("judge whether there are different behavior cyclones in same day for each lev")
for nl in range(0,len(lev),1):
    for nr1 in range(1,len(behv)-1,1):
        for nr2 in range(nr1+1,len(behv),1):
            ct1 = locals().get('ct_%d_%s'%(lev[nl],behv[nr1])) +\
                locals().get('ct_%d_%s'%(lev[nl],behv[nr2]))
            ct2=list(set(ct1))
            num=len(ct1)-len(ct2)
            num1=len(locals().get('ct_%d_%s'%(lev[nl],behv[nr1])))
            num2=len(locals().get('ct_%d_%s'%(lev[nl],behv[nr2])))
            if num1!=0 and num2!=0:
                print("repeated time in %dlev for %s and %s : %d, %.2f%%, %.2f%%"\
                        %(lev[nl],behv[nr1],behv[nr2],num,num*100/num1,num*100/num2))

print("")
print("combine 3levs ct for each behavior,delete repeated time")
print("total, unrepeated, repeated, repeated/total, unrepeated/alltime")
for nr in range(1,len(behv),1):
    locals()['ct_%s'%behv[nr]] = []
    for nl in range(0,len(lev),1):
        locals().get('ct_%s'%behv[nr]).extend(locals().get('ct_%d_%s'%(lev[nl],behv[nr])))
    
    locals()['ct2_%s'%behv[nr]] = list(set(locals().get('ct_%s'%behv[nr])))
    num1=len(locals().get('ct_%s'%behv[nr])) # total time
    num2=len(locals().get('ct2_%s'%behv[nr])) # unrepeated
    print("all time in three levs for %s: %d, %d, %d, %.2f%%, %.2f%%"\
            %(behv[nr],num1,num2,(num1-num2),(num1-num2)*100/num1,num2*100/alltime))

