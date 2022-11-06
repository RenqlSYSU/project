#!/usr/bin/env python
import sys, os, subprocess, linecache
import numpy as np
from datetime import datetime

def composite_time(filname,flats,flatn,flonl,flonr,nint):
    prefix = filname.split("_",1)[0].rsplit("/")[-1]
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    ctime=[]
    clat =[]
    clon =[]
    cinte=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            for nl in range(0,num,1):
                line = ff.readline()
                if prefix in ['ffadd','fftadd']:
                    data = list(map(float,line.strip().replace(" &","").split(" ")))
                else:
                    data = list(map(float,line.strip().split(" ")))
                if data[1]<=flonr and data[1] >= flonl and\
                data[2]<=flatn and data[2]>=flats :
                    ctime.append(datetime.strptime(str(int(data[0])),'%Y%m%d%H'))
                    clat.append(data[2])
                    clon.append(data[1])
                    cinte.append(data[nint])

        line = ff.readline()
    ff.close()
    return ctime,clat,clon,cinte

