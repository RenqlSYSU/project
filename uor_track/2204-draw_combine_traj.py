#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool
import sys, os, subprocess, linecache, gc
from datetime import datetime
from renql import cyc_filter, draw
import matplotlib.pyplot as plt
from matplotlib import colors,colorbar
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import shapely.geometry as sgeom
import cmaps

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

flats = 25  #int(sys.argv[2])
flatn = 45  #int(sys.argv[3])
flonl = 60  #int(sys.argv[4])
flonr = 110 #int(sys.argv[5])
prefix = "fftadd"
suffix = ['_0_2545-60105','_5_2545-60110_2_2545-6060','_5_2545-60110_2_4545-60110']
title= {'_0_2545-60105':'Local',
        '_5_2545-60110_2_2545-6060':'Western',
        '_5_2545-60110_2_4545-60110':'Northern'}
behv = {'_0_2545-60105':["PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS"],
        '_5_2545-60110_2_2545-6060':["PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ],
        '_5_2545-60110_2_4545-60110':["EPAS","ELYS","WPAS","WNTP","WNTL","WLYS"]}
lev  = [850,500,250]
colrs= ['b','g','r','c','m','y']

fileout="/home/users/qd201969/uor_track/mdata/"
path = '/home/users/qd201969/ERA5-1HR-lev/'
figdir = '/home/users/qd201969/uor_track/fig/'
datapath = "/work/scratch-pw2/renql/ERA5_hourly/wind10/ERA5_speed10_1hr_dec-jan"

lonl=0  #0  #
lonr=150#360#
lats=15  #
latn=70 #

def main_run():
    tlons = np.arange(0,150,1)
    draw_comb_traj(tlons)

def draw_comb_traj(tlons):
    lat_sp = 20
    lon_sp = 30 #60 #
    nrow = 3
    ncol = 3 
    bmlo = 0.45 #0.25 #
    title_font=14
    label_font=10
    
    ncmap = cmaps.precip2_17lev
    cnlvl1 = np.arange(0,8.5,0.5) 
    norm  = colors.BoundaryNorm(boundaries=cnlvl1, 
            ncolors=ncmap.N,extend="both")
    
    width = np.arange(0.1, 22.1, 2)
    cnlvl2 = np.arange(0, 1, 0.1)

    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    lat = ds.lat
    lon = ds.lon
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds

    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree()))
    
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            axe = ax[nr][nc]
            axe.set_xlim(lonl,lonr)
            axe.set_ylim(lats,latn)
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k')
                    , linewidth=0.8, zorder=1)
            axe.set_title('%s %dhPa'%(title[suffix[nc]],lev[nr]),
                fontdict=font, fontsize=title_font)
            
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                 linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            
            total = count_cyclone("%s%s_%d_1980-2020%s"%(path,prefix,lev[nr],suffix[nc]))
            for bhv,colr in zip(behv[suffix[nc]],colrs):
                filname = "%s%s_%d_1980-2020%s_%s"%(path,prefix,lev[nr],suffix[nc],bhv)
                if os.path.isfile("%s_ctraj"%filname):
                    data = np.loadtxt("%s_ctraj"%filname,delimiter=',')
                    if len(data)==0:
                        clon = []
                    else:
                        clon = list(data[:,0])
                        clat = list(data[:,1])
                        inte = list(data[:,2])
                        numb = list(data[:,3])
                else:
                    clon,clat,inte,numb = combine_traj(filname,tlons)

                if len(clon) > 2:
                    for i in range(len(clon)-1):
                        axe.add_geometries([sgeom.LineString([(clon[i],clat[i]),(clon[i+1],clat[i+1])])],
                            color=colr,linewidth=get_width(numb[i]/total,width,cnlvl2),crs=ccrs.PlateCarree())#
                            #get_color(inte[i],ncmap,cnlvl1)
                    axe.arrow(clon[-2],clat[-2],(clon[-1]-clon[-2]),(clat[-1]-clat[-2]),
                            color=colr,head_width=4,transform=ccrs.PlateCarree())
                    axe.plot(clon[0],clat[0],color=colr,label=bhv)
                
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
                axe.legend(loc='center',bbox_to_anchor=(0.5,-0.5),ncol=3) 
    #fig.legend(loc='center',bbox_to_anchor=(0.5,bmlo+0.005),ncol=6) 
    #position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
    #colourMap = plt.cm.ScalarMappable(cmap=ncmap,norm=norm)
    #fig.colorbar(colourMap,cax=position,extend='both',shrink=0.7,
    #        orientation='horizontal',
    #        label='Relative Vorticity (Hz) * 1e5')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig("%scombie_traj.png"%figdir,bbox_inches='tight',pad_inches=0.01)

def get_color(var,ncmap,cnlevels):
    for nll in range(0,len(cnlevels),1):
        if nll==0 and var<cnlevels[nll]:
            color0 = ncmap([nll])[0]
            break
        elif nll>0 and nll<=(len(cnlevels)-1) and var>=cnlevels[nll-1] and var<cnlevels[nll]:
            color0 = ncmap([nll])[0]
            break
        elif nll==(len(cnlevels)-1) and var>=cnlevels[nll]:
            color0 = ncmap([nll+1])[0]
    print(color0)
    return color0

def get_width(var,ncmap,cnlevels):
    for nll in range(0,len(cnlevels),1):
        if nll==0 and var<cnlevels[nll]:
            wh = ncmap[nll]
            break
        elif nll>0 and nll<=(len(cnlevels)-1) and var>=cnlevels[nll-1] and var<cnlevels[nll]:
            wh = ncmap[nll]
            break
        elif nll==(len(cnlevels)-1) and var>=cnlevels[nll]:
            wh = ncmap[nll+1]
    print(wh)
    return wh

def combine_traj(filname,tlons):
    clat = []
    clon = []
    inte = []
    numb = []
    ff = open("%s_ctraj"%filname, "w")
    for tlon in tlons:
        tlat,avginte = calc_average_traj(filname,tlon)
        # return two lists

        if len(tlat) > 5 :
            clon.append(tlon)
            numb.append(len(tlat))
            clat.append(avg_mvmaxmin(tlat))
            inte.append(avg_mvmaxmin(avginte))
            print("%.2f: %d cyclon, %.2f lat, %.2f vorticity"%(
                clon[-1],numb[-1],clat[-1],inte[-1]))
            ff.write("%f,%f,%f,%d \n"%(
                clon[-1],clat[-1],inte[-1],numb[-1]))
    ff.close() 
    return clon,clat,inte,numb

def count_cyclone(filname):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    total = a[1].strip().split(" ",1)[0]
    ff.close()
    return int(total)

def calc_average_traj(filname,tlon):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    total = a[1].strip().split(" ",1)[0]
    
    clat = [] 
    inte = []
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            clat1 = []
            inte1 = []
            for nl in range(0,num,1):
                line = ff.readline()
                data = list(map(float, line.strip().replace(" &","").split(" ")))
                if data[1]<(tlon+0.5) and data[1]>=(tlon-0.5):
                    clat1.append(data[2])
                    inte1.append(data[3])
            if len(clat1) > 0:
                clat.append(np.mean(clat1))
                inte.append(np.mean(inte1))
        line = ff.readline()
    ff.close()
    return clat,inte

def avg_mvmaxmin(var1d):
    if len(var1d) > 2:
        var1d.remove(max(var1d)), var1d.remove(min(var1d))
    return np.mean(var1d) 

if __name__=='__main__':
    main_run()
