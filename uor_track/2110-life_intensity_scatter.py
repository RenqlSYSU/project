#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-5960
then use box to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import sys, os, subprocess, linecache
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps

if len(sys.argv) < 2 :
    filname = "/home/users/qd201969/ERA5-1HR-lev/ff_250_1980-2020_2_3045-6060"
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

tid = []
life = []
inte = [] #  max intensity
tlat = []  # max intensity location
tlon = []  # max intensity location
line = ff.readline()
while line:
    term = line.strip().split(" ")
    if term[0] == "TRACK_ID":
        tid.append(term[-1])
        
        linenum = ff.readline()
        term1 =linenum.strip().split(" ")
        num = int(term1[-1])
        life.append(num/24.0)
        
        data=[]
        for nl in range(0,num,1):
            data.append(list(map(float, ff.readline().strip().split(" "))))

        data = np.array(data)
        inte.append(data[:,3].max())
        loc = np.argmax(data[:,3])
        tlat.append(data[loc,2])
        tlon.append(data[loc,1])

    line = ff.readline()

a = line4.strip().split(" ",1)
term = a[1].strip().split(" ",1)
print("total cyclone number in %s : %s" %(ff.name,term[0]))
ff.close()

print("tid life(days) intensity longitude latitude")
for ni in range(0,len(tid),1):
    print("%s "%tid[ni] + "%f "*4 %(life[ni],inte[ni],tlon[ni],tlat[ni]))

#===============================================
# draw figure 
#===================================================
lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30
nrow = 2
ncol = 1
bmlo = 0.4
BIGFONT=22
MIDFONT=14
SMFONT=10

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
lat = ds.lat
lon = ds.lon
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds

term = filname.split("/")
fig = plt.figure(figsize=(9,9),dpi=200)

axe = plt.subplot(nrow,ncol,1,projection=ccrs.PlateCarree())    #创建子图
axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
axe.set_title(term[-1],fontsize=MIDFONT)

axe.contour(ilon, ilat, phis, [1500,3000], 
      transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
axe.scatter(tlon,tlat,transform=ccrs.PlateCarree())
axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

ax = plt.subplot(nrow,ncol,2)    #创建子图
ax.scatter(life,inte)
ax.set_xlabel("lifetime (days)",fontsize=MIDFONT)
ax.set_ylabel("Max intensity ($10^{-5} s^{-1}$)",fontsize=MIDFONT)
ax.set_title(term[-1],fontsize=MIDFONT)

fig.savefig("/home/users/qd201969/uor_track/fig/life_inte_"+term[-1]+".png")

