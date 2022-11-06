#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use line to filter different behaviors cyclones

20210928
'''

import sys, os, subprocess, linecache

#def intersection_point(line1, line2):
#    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
#    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])
#
#    def det(a, b):
#        return a[0] * b[1] - a[1] * b[0]
#
#    div = det(xdiff, ydiff)
#    if div == 0:
#        raise Exception('lines do not intersect')
#
#    d = (det(*line1), det(*line2))
#    x = det(d, xdiff) / div
#    y = det(d, ydiff) / div
#    return x, y

def intersection_point_fixy(line1, y):
    xdiff = line1[0][0] - line1[1][0]
    ydiff = line1[0][1] - line1[1][1]
    if xdiff == 0 :
        x = line1[0][0]
    else:
        slope = ydiff / xdiff
        intercept = line1[0][1] - slope*line1[0][0]
        x = (y-intercept)/slope
    return x

def intersection_point_fixx(line1, x):
    xdiff = line1[0][0] - line1[1][0]
    ydiff = line1[0][1] - line1[1][1]
    slope = ydiff / xdiff
    intercept = line1[0][1] - slope*line1[0][0]
    y = slope*x + intercept
    return y

#if sys.argv[1] == "-1" :
if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 0 #int(sys.argv[2])
    flatn = 90 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 60 #int(sys.argv[5])
    time = 24 # threshold, hour
    filname = "/home/users/qd201969/ERA5-1HR-lev/ff_850_1980-2020"
else:
    filname = sys.argv[1]
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    option= int(sys.argv[6]) 
    time = int(sys.argv[7]) 
suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
print(filname+"_"+suffix+'\n')

#===============================================
# read ff_trs_pos file and write output file
#===================================================
ff = open(filname,"r") 
outfile = open(filname+"_"+suffix,"w")
line1 = ff.readline()
line2 = ff.readline()
line3 = ff.readline()
line4 = ff.readline()
outfile.write(line1)
outfile.write(line2)
outfile.write(line3)
outfile.write(line4)

tid=[]
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
            data.append(list(map(float, line.strip().split(" "))))
        
        if num >= time :
            if flonl == flonr :
                for nl in range(0,num-1,1):
                    if data[nl][1] <= flonl and data[nl+1][1] >= flonl and data[nl+1][1] <= 180 :
                        point = intersection_point_fixx([data[nl][1:3], data[nl+1][1:3]], flonl)
                        if point <= flatn and point >= flats :
                            #print(lineid.strip()+" " +linenum.strip()+" pass "+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr))
                            tid.append(term[-1])
                            outfile.write(lineid)
                            outfile.write(linenum)
                            #outfile.write(str(data).strip('[').strip(']').replace('], [','\n').replace(',','')+'\n')
                            for nll in value:
                                outfile.write(nll)
                            break 

            if flats == flatn :
                for nl in range(0,num-1,1):
                    if (data[nl][2]-flats)*(data[nl+1][2]-flats) <= 0 :
                        point = intersection_point_fixy([data[nl][1:3], data[nl+1][1:3]], flats)
                        if point <= flonr and point >= flonl :
                            #print(lineid.strip()+" " +linenum.strip()+" pass "+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr))
                            tid.append(term[-1])
                            outfile.write(lineid)
                            outfile.write(linenum)
                            #outfile.write(str(data).strip('[').strip(']').replace('], [','\n').replace(',','')+'\n')
                            for nll in value:
                                outfile.write(nll)
                            break 

    line = ff.readline()

outfile.seek(0,0) # Go back to the beginning of the file
outfile.write(line1)
outfile.write(line2)
outfile.write(line3)
outfile.write("TRACK_NUM %9d ADD_FLD    0   0 &" %len(tid))

a = line4.strip().split(" ",1)
term = a[1].strip().split(" ",1)
print("total cyclone number in %s : %s" %(ff.name,term[0]))
print("filter cyclone number in %s : %d" %(outfile.name,len(tid)))
ff.close()
outfile.close()
#print(tid)

