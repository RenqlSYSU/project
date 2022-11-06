#!/usr/bin/env python
'''
read yearly cyclone in ff_trs_pos 
output cyclone exist in 12 to 11 month
to avoid overlap time periods
then transfer timestep to real time

20220120
renql
'''

import numpy as np
import sys, os, subprocess

def month_filter(filname):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    a = line3.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    outfile = open("%s-1211"%(filname),"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    
    tid=[]
    line = ff.readline()
    while line:
        term3 = line.strip().split(" ")
        if term3[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])

            ct1=[]
            value=[]
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
                data = list(map(float, line.strip().replace(" &","").split(" ")))
                #data = list(map(float,line.strip().split(" ")))
                ct1.append(int(data[0]))

            if sum(i<=9504 and i>=745 for i in ct1)/len(ct1) >= 0.5 : 
                tid.append(term3[2])
                outfile.write(lineid)
                outfile.write(linenum)
                for nll in value:
                    outfile.write(nll)
        
        line = ff.readline()
    
    ff.close()
    outfile.seek(0,0) # Go back to the beginning of the file
    outfile.write(line1)
    outfile.write(line2)
    outfile.write("TRACK_NUM %9d %s" %(len(tid),term[-1]))
    print("%s : %d" %(outfile.name,len(tid)))
    outfile.close()
    
    return len(tid)

path = '/home/users/qd201969/ERA5-1HR-lev/'
lev  = [850,500,250]
prefix = 'ff'

for year in range(1980,2021):
    for nl in lev:
        filname = "%sERA5_VOR%d_1hr_%d_DET/%s_trs_pos.addZ_addwind10m_precip_maxpre"\
                %(path,nl,year,prefix)
        if not os.path.exists('%s-1211'%filname):
            month_filter(filname)

com = "sh ~/uor_track/control_era5_1hr_track.sh %s -1 0 addZ_addwind10m_precip_maxpre-1211"%(prefix)
ret=subprocess.Popen(com,shell=True)
ret.wait()

com = "sh ~/uor_track/control_era5_1hr_track.sh %st 1 0 addZ_addwind10m_precip_maxpre-1211"%(prefix)
ret=subprocess.Popen(com,shell=True)
ret.wait()

