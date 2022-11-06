#!/usr/bin/env python
import sys, os, subprocess, linecache
import numpy as np

def filt_close_cyclone(filname,perc,flats=25,flatn=45,flonl=57,flonr=110):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term = a[1].strip().split(" ",1)
    print("total cyclone number in %s : %s" %(ff.name,term[0]))
    
    outfile = open("%s_%dclose"%(filname,perc),"w")
    outfile.write(line1)
    outfile.write(line2)
    outfile.write(line3)
    outfile.write(line4)

    tid = []
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            lineid = line
            linenum = ff.readline()
            term1 = linenum.strip().split(" ")
            num = int(term1[-1])
            
            data = []
            value=[]
            clse_time=0
            for nl in range(0,num,1):
                line = ff.readline()
                value.append(line)
                term = list(map(float,line.strip().replace(" &","").split(" ")))
                data.append(term)

            signal = 0
            data = np.array(data)
            clse1 = 100-sum(i==1e+25 for i in data[:,4])*100/len(data)
            if clse1 >= perc:
                signal = 1
                print("%d close cyclone: %d-%d lon:%.2f-%.2f lat:%.2f-%.2f"%(
                    perc, data[0,0],data[-1,0],data[0,1],data[-1,1],data[0,2],data[-1,2]))

            if signal == 1:
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
    outfile.write("TRACK_NUM %9d ADD_FLD    0   0 &" %len(tid))
    print("filter cyclone number in %s : %d" %(outfile.name,len(tid)))
    outfile.close()

    return len(tid)

