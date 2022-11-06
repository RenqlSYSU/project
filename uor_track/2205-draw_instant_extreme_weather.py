#!/usr/bin/env python
'''
fix a time and then draw the instant geopotential (contour) from 
/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/ERA5_NH_z_1989.nc,

masked extreme 10m wind (shaded) from 
~/ERA5-1HR-lev/ERA5_VOR850_1hr_1995_DET/ERA5_VOR850_1hr_1995_DET_T63filt.nc

and identified feature points from 
~/ERA5-1HR-lev/ff_250_500_yes_850_yes_match
ff_500_250_yes_850_yes_match
ff_850_500_yes_250_yes_match

20220515
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
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

dtime = pd.date_range(start='1995-01-04 00',periods=60, freq='12H',closed=None)
#dtime = pd.date_range(start='1995-01-01 00',end='1995-01-15 00', freq='6H',closed=None)
lev = [850,500,250]
cnlvl2 = [30,50,100]
varname = 'z'
path = '/home/users/qd201969/ERA5-1HR-lev/'
datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily"#t/ERA5_NH_t_1989.nc
figdir = "/home/users/qd201969/uor_track/fig/"
filname=['fftadd_850_1980-2020','fftadd_500_1980-2020','fftadd_250_1980-2020']
#filname=['ff_850_500_yes_250_yes_match',
#    'ff_500_250_yes_850_yes_match','ff_250_500_yes_850_yes_match']
radiu = 6

def main_run():
    draw_instant()

def draw_instant():
    f  = xr.open_dataset("%s/%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,dtime[0].year))
    lat = f['latitude'].data
    lon = f['longitude'].data
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
    phis = phis/9.8 # transfer from m2/s2 to m
    del ds
    gc.collect()
    
    #ds = xr.open_dataset("/home/users/qd201969/uor_track/mdata/max10mwind_99.0threshold_month.nc")
    ds = xr.open_dataset("/home/users/qd201969/uor_track/mdata/maxprecip1h_99threshold_month.nc")
    ilon1 = ds.lon
    ilat1 = ds.lat
    thre = ds['threshold'].data

    ncmap = colors.ListedColormap(["blue","red","green"])
    cnlevels = [0.5,2]
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=3,extend="both")

    for nt in range(len(dtime)):
        fig = plt.figure(figsize=(12,12),dpi=100)
        ax = fig.subplots(nrow,ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
        #fvor = xr.open_dataset("/work/scratch-pw2/renql/ERA5_hourly/wind10/ERA5_speed10_1hr_dec-jan%d.nc"%(dtime[nt].year))
        #var1 = fvor['var1'].sel(time=dtime[nt]).load()
        fvor = xr.open_dataset("/gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/precip/ERA5_precip_1hr_dec-jan%d.nc"%(dtime[nt].year))
        var1 = fvor['tp'].sel(time=dtime[nt],longitude=ilon1,latitude=ilat1).load()
        #var1 = xr.where( var1>thre[dtime[0].month-1,:,:], 1, 0 )
        var1 = xr.where( var1>0, 1, 0 )
        print('extreme weather point: %d'%(np.sum(var1.data==1)))

        for nl in range(len(lev)):
            var = f[varname].sel(time=dtime[nt],level=lev[nl],longitude=ilon,latitude=ilat)
            var.data = var.data/9.8
            plat, plon = read_point_fixtime("%s/%s"%(path,filname[nl]), 
                dtime[nt].strftime('%Y%m%d%H'),lonl,lonr,lats,latn)
            
            mask = np.zeros(var1.shape,dtype=float)
            for clo,cla in zip(plon,plat):
                mask = np.where((np.square(var1.lat-cla)+np.square(var1.lon-clo))<=(radiu*radiu),1,mask)
            term = np.ma.array(var1.data, mask=(mask==0))
            
            axe = ax[nl]
            axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            axe.set_title("%s %dhPa (%d)"%(dtime[nt].strftime('%Y-%m-%d-%H:00'), lev[nl], len(plat)),fontsize=title_font)

            shad = axe.contourf(var1.lon, var1.lat, term, cnlevels,
                    transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
            
            cont = axe.contour(ilon, ilat, var, np.arange(1000,15000,cnlvl2[nl]), 
                    transform=ccrs.PlateCarree(), colors='gray', linewidths=1.5)
            
            pint = axe.scatter(plon,plat,10.0**2,color='k', marker='o', transform=ccrs.PlateCarree())

            topo = axe.contour(ilon, ilat, phis, [1500,3000],
                    transform=ccrs.PlateCarree(),colors='black',linewidths=1.2)

            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

        plt.tight_layout(rect=(0,bmlo,1,1))
        plt.savefig(figdir+"maxprecip_%s.png"%(dtime[nt].strftime('%Y%m%d%H')), bbox_inches='tight',pad_inches=0.01)

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
                #data = list(map(float,ff.readline().strip().split(" ")))
                data = list(map(float, ff.readline().strip().replace(" &","").split(" ")))
                if str(int(data[0])) == fixtime and \
                data[1]<=flonr and data[1] >= flonl and data[2]<=flatn and data[2]>=flats :
                    plat.append(data[2])
                    plon.append(data[1])
        line = ff.readline()
    ff.close()
    print("%s total feature point in %s : %d"%(filname,fixtime,len(plat)))
    return plat, plon 

def create_gif(figname):
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

if __name__=='__main__':
    main_run()
