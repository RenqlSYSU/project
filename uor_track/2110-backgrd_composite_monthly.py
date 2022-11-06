#!/usr/bin/env python
'''
1. when cyclone located in the defined region, read datetime
2. then convert hourly date to 6 hourly date,
3. composite corresponding u,v,t,z and their variances 
4. store the data

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
from renql import cyc_filter

def mean_vari_dynamic(mean0,vari0,numb0,samp):
    print("numb0 %d, mean0[1,12,12] %.2f, vari0[1,12,12] %.2f"%(numb0,mean0[1,12,12],vari0[1,12,12]))
    numb1 = numb0 + len(samp)
    mean1 = mean0 + (np.sum(samp,axis=0)-len(samp)*mean0)/numb1
    vari1 = (numb0*(vari0+np.power(mean1-mean0,2))+np.sum(np.power(np.subtract(samp,mean1),2),axis=0))/numb1
    return mean1,vari1,numb1

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
    inty=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data = []
            for nl in range(0,num,1):
                line = ff.readline()
                data.append(list(map(float,line.strip().split(" "))))
                '''
                data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    a = str(int(data[0]))
                '''
            for nl in range(0,num-1,1):
                if data[nl][1] <= flonl and data[nl+1][1] > flonl and data[nl+1][1] <= 180 :
                    point = cyc_filter.intersection_point_fixx([data[nl][1:3], data[nl+1][1:3]], flonl)
                    if point <= flatn and point >= flats :
                        a = str(int(data[nl+1][0]))
                        '''
                        if int(a[-2:]) in range(0,6,1):
                            b='00'
                        elif int(a[-2:]) in range(6,12,1):
                            b='06'
                        elif int(a[-2:]) in range(12,18,1):
                            b='12'
                        elif int(a[-2:]) in range(18,24,1):
                            b='18'
                        '''
                        a = a[:-2]+'00'
                        ct1.append(datetime.strptime(a,'%Y%m%d%H'))
                        inty.append(data[nl+1][3])
        line = ff.readline()
    ff.close()
    
    print("total, unrepeated, repeated, repeated/total, unrepeated/alltime")
    if len(ct1) == 0:
        ct2 = ct1
        print("all time when cyclone in %d-%dE,%d-%dN : %d"\
                %(flonl,flonr,flats,flatn,len(ct1)))
    else:
        #ct1=list(np.array(ct1)[inty>np.percentile(np.array(inty),95)]) # top 5%
        #ct1=list(np.array(ct1)[inty<np.percentile(np.array(inty),5)]) # low 5%
        ct2=list(set(ct1))
        print("all time when cyclone in %d-%dE,%d-%dN : %d, %d, %d, %.2f%%, %.2f%%"\
                %(flonl,flonr,flats,flatn,len(ct1),len(ct2),(len(ct1)-len(ct2)),\
                (len(ct1)-len(ct2))*100/len(ct1),len(ct2)*100/alltime))
    print("")
    return ct1

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 25 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 105 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "fft"
    suffix = '5_2545-60110_2_2545-6060'
    behv = ["ALL" ,"PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
    time = int(sys.argv[8])

fileout="/home/users/qd201969/uor_track/mdata/comp_6h_season_daily00_"
levc = [850,500,250,200]
lev  = [850,500,250]
path = '/home/users/qd201969/ERA5-1HR-lev/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/"#t/ERA5_NH_t_1989.nc

ftime  = pd.date_range(start='1979-12-01 00',end='2021-01-31 23', freq='1D',closed=None)
alltime= len(ftime)

for nl in range(0,len(lev),1):
    for nr in range(1,len(behv),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix+"_"+behv[nr]
        locals()['ct_%d_%s'%(lev[nl],behv[nr])]=composite_time(filname,flats,flatn,flonl,flonr,alltime)

months = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
season = ["DJF","MAM","JJA","SON"]

lonl=0  #0  #
lonr=150#360#
lats=0  #
latn=90 #
ds  = xr.open_dataset(datapath+"t/ERA5_NH_t_1989.nc")
lat = ds.latitude
lon = ds.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
for varname in ['u','v','z']: #'t',
    var  = np.zeros( [len(lev),len(behv),len(months),len(levc),len(ilat),len(ilon)],dtype=float )  
    vari = np.zeros( [len(lev),len(behv),len(months),len(levc),len(ilat),len(ilon)],dtype=float )  
    numb = np.zeros( [len(lev),len(behv),len(months)],dtype=int ) 
    for nl in range(0,len(lev),1):
        for nr in range(0,len(behv),1):
            print("handle %d %s"%(lev[nl],behv[nr]))
            if nr == 0:
                for nm in range(0,len(months),1):
                    term_mean = np.zeros( [len(levc),len(ilat),len(ilon)],dtype=float ) 
                    term_vari = np.zeros( [len(levc),len(ilat),len(ilon)],dtype=float ) 
                    allt = 0
                    for year in range(1979,2021,1):
                        ds   = xr.open_dataset("%s%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,year))
                        term = ds[varname].sel(time=np.array([ds.time.dt.month.isin(months[nm]),ds.time.dt.hour.isin(0)]).all(axis=0),
                                level=levc,longitude=ilon,latitude=ilat)
                        term_mean, term_vari, allt = mean_vari_dynamic(term_mean,term_vari,allt,term)
                        #allt = allt + len(ds['time'].sel(time=ds.time.dt.month.isin(months[nm])))
                    var[nl,nr,nm,:,:,:] = term_mean
                    vari[nl,nr,nm,:,:,:]= term_vari
                    numb[nl,nr,nm]=allt
                    del term, allt, term_mean, term_vari
            else:
                ct = locals().get('ct_%d_%s'%(lev[nl],behv[nr]))
                if len(ct) == 0:
                    numb[nl,nr,:]= 0
                    var[nl,nr,:,:,:,:] = 0 
                    vari[nl,nr,:,:,:,:] = 0 
                    continue
                ctda = xr.DataArray(ct, coords=[ct], dims=["time"])
                for nm in range(0,len(months),1):
                    term_mean = np.zeros( [len(levc),len(ilat),len(ilon)],dtype=float ) 
                    term_vari = np.zeros( [len(levc),len(ilat),len(ilon)],dtype=float ) 
                    allt = 0
                    ctda0 = ctda.sel(time=ctda.time.dt.month.isin(months[nm]))
                    if len(ctda0) == 0:
                        numb[nl,nr,nm]=0
                        var[nl,nr,nm,:,:,:] = 0 
                        vari[nl,nr,nm,:,:,:] = 0 
                        continue
                    for year in range(1979,2021,1):
                        ctda1 = ctda0.sel(time=ctda0.time.dt.year.isin(year))
                        #print("%d numb %d"%(year,len(ctda1)))
                        if len(ctda1) == 0:
                            continue
                        ds  = xr.open_dataset("%s%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,year))
                        term = ds[varname].sel(time=ctda1,level=levc,longitude=ilon,latitude=ilat)
                        term_mean, term_vari, allt = mean_vari_dynamic(term_mean,term_vari,allt,term)
                    var[nl,nr,nm,:,:,:] = term_mean
                    vari[nl,nr,nm,:,:,:]= term_vari
                    numb[nl,nr,nm]=allt
                del ct, ctda, term, ctda0, ctda1, allt, term_mean, term_vari
            
    ds = xr.Dataset(
            {
                "numb": (["lev3", "behv","month"], numb),
                "var" : (["lev3", "behv","month", "level", "lat","lon"], var),
                "vari": (["lev3", "behv","month", "level", "lat","lon"], vari),
                },
            coords={
                "lev3" : (["lev3"],lev),
                "behv" : (["behv"],behv),
                "month": (["month"],season),
                "level": (["level"],levc),
                "lat"  : (["lat"],ilat.data),
                "lon"  : (["lon"],ilon.data),
                },
            )
    ds.attrs["description"] = "composite background for diff behavior of 3levs cyclones"
    ds.to_netcdf(fileout+varname+".nc","w")
    del ds,var,numb,vari
        
