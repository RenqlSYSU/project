#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use box to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import sys, os, subprocess, linecache

if len(sys.argv) < 2 :
    filname = "/home/users/qd201969/ERA5-1HR-lev/ff_850_1980-2020_2_3045-6060"
else:
    filname = int(sys.argv[1])

#===============================================
# read ff_trs_pos file and write output file
#===================================================
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
#        print(lineid.strip() + " " + linenum.strip())
        
        data=[]
        for nl in range(0,num,1):
            data.append(list(map(float, ff.readline().strip().split(" "))))

        imax = 0
        imin = 0
        inl = []
        for nl in range(1,num-1,1):
            tend1 = data[nl][3] - data[nl-1][3]
            tend2 = data[nl+1][3] - data[nl][3]
            if tend1>=0 and tend2<=0:
                imax += 1
                inl.append(nl)
            if tend1<=0 and tend2>=0:
                imin += 1
                inl.append(nl)
        
        if len(inl)>=2:
            print(lineid.strip()+" "+linenum.strip()+" has imax: %d, imin: %d" %(imax,imin)) 
            for nl in inl:
                print(str(nl)+" "+str(data[nl]).strip('[').strip(']').replace(',',''))
            print('\n')
            tid.append(term[-1])

    line = ff.readline()

a = line4.strip().split(" ",1)
term = a[1].strip().split(" ",1)
print("total cyclone number in %s : %s" %(ff.name,term[0]))
print("%d cyclone have multiple extreme intensity" %len(tid))
ff.close()


