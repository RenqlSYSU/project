#!/usr/bin/env python
'''
fix a time and then draw the instant geopotential (contour) from 
/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/ERA5_NH_z_1989.nc,

spatial filtered relative vorticity (shaded) from 
~/ERA5-1HR-lev/ERA5_VOR850_1hr_1995_DET/ERA5_VOR850_1hr_1995_DET_T63filt.nc

and identified feature points from 
~/ERA5-1HR-lev/ERA5_VOR850_1hr_1995_DET/fft_trs_pos

Loop through the height (850, 500, 250)

20211116
'''
import sys
import subprocess
import xarray as xr
import numpy as np
import pandas as pd
from datetime import datetime
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from PIL import Image, ImageDraw, ImageSequence

def calc_frames(new_time):
    old_time = datetime(new_time.year-1, 11, 30, 23)
    days = (new_time - old_time).days
    sec = (new_time - old_time).seconds
    hours = days * 24 + sec/3600
    return int(hours)

def read_point_fixtime(filname,fixtime,flonl,flonr,flats,flatn):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    
    plat = []
    plon = []
    line = ff.readline()
    while line:
        if line.strip().split(" ")[0] == "TRACK_ID":
            num = int(ff.readline().strip().split(" ")[-1])
            for nl in range(0,num,1):
                data = list(map(float,ff.readline().strip().split(" ")))
                if str(int(data[0])) == fixtime and \
                data[1]<=flonr and data[1] >= flonl and data[2]<=flatn and data[2]>=flats :
                    plat.append(data[2])
                    plon.append(data[1])
        line = ff.readline()
    ff.close()
    print("%s total feature point in %s : %d"%(filname,fixtime,len(plat)))
    return plat, plon 

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30

nrow = 3
ncol = 1
bmlo = 0.1
title_font=18
label_font=14

dtime = pd.date_range(start='1995-01-01 00',periods=60, freq='6H',closed=None)
#dtime = pd.date_range(start='1995-01-01 00',end='1995-01-15 00', freq='6H',closed=None)
create_gif = True #False#
nfilt="T63"
lev = [850,500,250]
cnlvl =[[-8 ,1 ]]
cnlvl2 = [30,50,100]
varname = 'z'
path = '/home/users/qd201969/ERA5-1HR-lev/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/"#t/ERA5_NH_t_1989.nc
figdir = "/home/users/qd201969/uor_track/fig/"

f  = xr.open_dataset("%sERA5_subdaily/%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,dtime[0].year))
lat = f['latitude'].data
lon = f['longitude'].data
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

nl = 0
fcolors = cmaps.BlueDarkRed18
cnlevels = np.arange(cnlvl[nl][0], cnlvl[nl][0]+cnlvl[nl][1]*(fcolors.N-1), cnlvl[nl][1])
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

params = {'legend.fontsize': label_font,
          'axes.labelsize': label_font,
          'axes.titlesize':label_font,
          'xtick.labelsize':label_font,
          'ytick.labelsize':label_font}
plt.rcParams.update(params)

for nt in range(len(dtime)):
    fig = plt.figure(figsize=(12,12),dpi=100)
    ax = fig.subplots(nrow,ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
    for nl in range(len(lev)):
        var = f[varname].sel(time=dtime[nt],level=lev[nl],longitude=ilon,latitude=ilat)
        var.data = var.data/9.8

        path2 = "%sERA5_VOR%d_1hr_%d_DET/"%(path,lev[nl],dtime[nt].year)
        plat, plon = read_point_fixtime(path2+"fft_trs_pos",dtime[nt].strftime('%Y%m%d%H'),lonl,lonr,lats,latn)
        
        fvor = xr.open_dataset("%sERA5_VOR%d_1hr_%d_DET_%sfilt.nc"%(path2,lev[nl],dtime[nt].year,nfilt))
        var1 = fvor['var'].sel(time=calc_frames(dtime[nt]),level = 1,lon=ilon,lat=ilat,method="nearest").load()
        #fvor = xr.open_dataset("%sERA5_VOR_1h_dec_jan/ERA5_VOR%d_1hr_dec-jan%d_DET.nc"%(datapath,lev[nl],dtime[nt].year))
        #var1 = fvor['var138'].sel(time=dtime[nt],lev=float(lev[nl]*100),lat=ilat,lon=ilon,method="nearest").load()
        var1.values = var1.values*1e5

        axe = ax[nl]
        axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
        axe.set_title("%s %dhPa (%d)"%(dtime[nt].strftime('%Y-%m-%d-%H:00'), lev[nl], len(plat)),fontsize=title_font)

        shad = axe.contourf(ilon, ilat, var1, cnlevels,
                transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        
        cont = axe.contour(ilon, ilat, var, np.arange(1000,15000,cnlvl2[nl]), 
                transform=ccrs.PlateCarree(), colors='gray', linewidths=1.5)
        
        #pint = axe.plot(plon,plat,color='darkviolet', marker='o', markersize=12, transform=ccrs.PlateCarree())
        pint = axe.scatter(plon,plat,10.0**2,color='k', marker='o', transform=ccrs.PlateCarree())

        topo = axe.contour(ilon, ilat, phis, [1500,3000],
                transform=ccrs.PlateCarree(),colors='black',linewidths=1.2)

        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
        axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    position = fig.add_axes([0.85, bmlo+0.1, 0.015, 0.7]) #left, bottom, width, height
    cb = plt.colorbar(shad, cax=position ,orientation='vertical')#, shrink=.9)
    cb.set_label(label='T5~63 Relative Vort (1e5)', size=label_font) #, weight='bold'

    plt.tight_layout(rect=(0,bmlo,1,1))
    plt.savefig(figdir+"filt_vor_%s.png"%(dtime[nt].strftime('%Y%m%d%H')), bbox_inches='tight',pad_inches=0.01)

if create_gif == True:
    figname = figdir+"filt_vor_*.png"
    fn_stream = subprocess.check_output("ls "+figname, shell=True).decode('utf-8')
    fn_list   = fn_stream.split()
    print(fn_list[0])
    print('filenumber : '+str(len(fn_list)))
    gif_name = figname.rsplit("_",1)[0]+".gif" 

    frames = []
    for itm in fn_list:
        frame = Image.open(itm)
        frames.append(frame)

    frames[0].save(gif_name, save_all=True, append_images=frames[1:],\
                duration = 1000, loop=0, disposal=1)
    subprocess.run('rm -f %s'%(figname),shell=True)

