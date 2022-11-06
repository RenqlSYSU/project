#!/usr/bin/env python
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
#matplotlib.use('Agg')

lonl=15  #0  #
lonr=145#360#
lats=15 #20 #
latn=70 #90 #
lat_sp = 20
lon_sp = 30 #60 #
nrow = 4 #6 #
ncol = 3 #2 #
bmlo = 0.35#0.4 
title_font=14
label_font=10

plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        'size': title_font,
        }

# level == 0 #small range
cnlvl=[[0    ,320  ,20  ], # 0Feature Density
     [0    ,3.2  ,0.2 ], # 1Genesis Density
     [0    ,3.2  ,0.2 ], # 2Lysis Density
     [-1   ,1    ,1   ], # 3Mean Area
     [-0.8 ,0.8  ,0.1 ], # 4Mean Growth/Decay Rate
     #[0    ,1.6  ,0.1 ], # 4Mean Growth/Decay Rate
     [-1   ,1    ,1   ], # 5Mean Anisotropy
     [0    ,8    ,0.5 ], # 6Mean Lifetime
     [0    ,80   ,5   ], # 7Mean Speed
     [0    ,12.8 ,0.8 ], # 8Mean Intensity
     #[0    ,8  ,0.5 ], # 9Mean Tendency
     [-1.6 ,1.6  ,0.2 ], # 9Mean Tendency
     [-1   ,1    ,1   ], # 10Spare1
     [-1   ,1    ,1   ], # 11Spare2
     [0    ,80   ,5   ], # 12Std of Speed
     [0    ,3.2  ,0.2 ], # 13Std of Intensity
     [0    ,16   ,1   ], # 14Track Density
     [-1   ,1    ,1   ], # 15X-component of Mean Orientation Vector
     [-40  ,40   ,5   ], # 16X-component of Mean Velocity
     [-1   ,1    ,1   ], # 17Y-component of Mean Orientation Vector
     [-40  ,40   ,5  ]] # 18Y-component of Mean Velocity

draw_var = ["fden","gden","lden","marea","mgdr","",
            "mlif","msp" ,"mstr","mten" ,""    ,"",
            ""    ,""    ,"tden",""    ] # 7 variables
draw=[7,8,]
#draw=[4,9,]
#draw=[6,7,8,4,9]
#draw=[1,2,14]
#draw=[14]
#draw=[1,2,14,8,9,6]
lev = [850,500,250]

if len(sys.argv) < 2 :
    prefix = "ff"
    suffix = ''
    dbox = 0 
    level = 2
else:
    prefix = sys.argv[1]  #'ff_250_500_no'
    suffix = sys.argv[2]  #'ff_250_500_no'
    level = int(sys.argv[3])
    dbox  = int(sys.argv[4])

if dbox >= 1 :
    flats = int(sys.argv[6])
    flatn = int(sys.argv[7])
    flonl = int(sys.argv[8])
    flonr = int(sys.argv[9])

if level == 1: #middle range 
    cnlvl[ 1][:]=[0    ,4.8  ,0.3 ]
    cnlvl[ 2][:]=[0    ,4.8  ,0.3 ]
    cnlvl[14][:]=[0    ,24   ,1.5 ]
if level == 2: # large range,use total level bar 
    cnlvl[ 0][:]=[0    ,1600 ,100 ]
    cnlvl[ 1][:]=[0    ,8    ,0.5 ]
    cnlvl[ 2][:]=[0    ,8    ,0.5 ]
    cnlvl[14][:]=[0    ,32   ,2   ]
titls=['DJF','MAM','JJA','SON']

files = '/home/users/qd201969/ERA5-1HR-lev/statistic/ff_250_1980-2020_stat.nc'
f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset('/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
# increased performance by loading data into memory first, e.g., with load()
uwnd = da.groupby(da.time.dt.month).mean('time')
print(uwnd)
del ds, da

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

figdir = "/home/users/qd201969/uor_track/fig/stat_season"+suffix
for nv in range(0,len(draw),1):#,len(f),1):
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
    
    cnlevels = np.arange(cnlvl[draw[nv]][0], cnlvl[draw[nv]][1]+cnlvl[draw[nv]][2], cnlvl[draw[nv]][2])
    if cnlvl[draw[nv]][0] < 0 :
        fcolors = cmaps.BlueDarkRed18
    else:
        fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
    
    for nl in range(0,len(lev),1):
        files = '/home/users/qd201969/ERA5-1HR-lev/statistic/%s_%d_1980-2020%s_stat.nc'%(prefix,lev[nl],suffix)
        print(files)
        f = xr.open_dataset(files)
        var = f[draw_var[draw[nv]]].sel(long=ilon,lat=ilat).load()
        if draw[nv] in [4,9]:
            var.data=var.data*24
        if draw[nv] > 2 and draw[nv] != 14:
            tden = f['tden'].sel(long=ilon,lat=ilat).load()
            mask = tden < 1.0
            var.values=np.ma.array(var.values,mask=mask)
    
        for nm in range(0,nrow,1):
            if nm == 0:
                var1 = (var[0,:,:]+var[1,:,:]+var[11,:,:])/3.0
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                var1 = np.mean(var[(3*nm-1):(3*nm+2),:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
            print('%s: min:%f ; max:%f'%(var.long_name,
                np.nanmin(var.data),np.nanmax(var.data)))
            axe = ax[nm][nl]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title("%dhPa %s"%(lev[nl],titls[nm]),fontdict=font)

            cont = axe.contourf(ilon, ilat, var1, cnlevels, 
                         transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                         transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            if dbox >= 1 :
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                         linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], 
                         transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2)
            if nl == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    #fig.subplots_adjust(left=0.2,bottom=bmlo) # wspace control horizontal space
    if (nrow/ncol) >= 2.0: 
        position = fig.add_axes([0.99, bmlo+0.05, 0.01, 0.6]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
        cb.set_label(label=var.long_name, size=title_font) #, weight='bold'
    else:
        position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
        plt.figtext(0.02,bmlo-0.005, var.long_name,fontsize=title_font,
                horizontalalignment='left',verticalalignment='bottom')
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir+draw_var[draw[nv]]+".png", bbox_inches='tight',pad_inches=0.01)


