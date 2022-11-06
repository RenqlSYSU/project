#!/usr/bin/env python
import sys
import subprocess, os
import xarray as xr
import numpy as np
import pandas as pd 
import gc #garbage collector
from datetime import datetime, timedelta
#matplotlib.use('Agg')

path = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev'
outdir = '/home/ys17-23/Extension2/renql/project/uor_track/mdata'
figdir = '/home/ys17-23/Extension2/renql/project/uor_track/fig'

def main_run():
    fn_stream = subprocess.check_output('ls %s/fftadd_match_*_6dist'%(
        outdir), shell=True).decode('utf-8')
    fn_list = fn_stream.split()
    for filname in fn_list:
        fftadd2ff(filname)

def fftadd2ff(filname):
    outfilename = filname.replace('fftadd','ff',1)
    if os.path.exists(outfilename):
        print('%s exists'%outfilename)
        return
    else:
        print('handle %s'%outfilename)
    
    ff = open(filname,"r") 
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

    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            outfile.write(lineid.strip().rsplit(" ",2)[0]+'\n')
            outfile.write(linenum)

            start = datetime.strptime(term[-1],'%Y%m%d%H')
            if start.month == 12:
                ftime = datetime(start.year, 12, 1, 00)
            else:
                ftime = datetime(start.year-1, 12, 1, 00)
            for nl in range(0,num,1):
                data = list(ff.readline().strip().replace(" &","").split(" "))
                data[0] = str(int((datetime.strptime(data[0],'%Y%m%d%H')-ftime
                    ).total_seconds()/60.0/60.0))
                outfile.write(' '.join(data[0:4])+'\n')

        line = ff.readline()
    ff.close()
    outfile.close()

if __name__=='__main__':
    main_run()

