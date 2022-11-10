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

lonl=0  #0  #
lonr=150#360#
lats=15 #10 #
latn=70 #90 #
lat_sp = 20 #30
lon_sp = 30 #60
nrow = 3
ncol = 3
bmlo = 0.45 #0.55 #
title_font=14
label_font=10
xbar=[0.05,0.37,0.69]
numod= [chr(i) for i in range(97,115)]

plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        'size': title_font,
        }

# level == 0 #small range
cnlvl=[[0    ,20  ], # 0Feature Density
       [0    ,0.2 ], # 1Genesis Density
       [0    ,0.2 ], # 2Lysis Density
       [-1   ,1   ], # 3Mean Area
       [-0.8 ,0.1 ], # 4Mean Growth/Decay Rate
       [-1   ,1   ], # 5Mean Anisotropy
       [2.2  ,0.2 ], # 6Mean Lifetime
       [20   ,3   ], # 7Mean Speed
       [2.5  ,0.5 ], # 8Mean Intensity
       [-1.6 ,0.2 ], # 9Mean Tendency
       [-1   ,1   ], # 10Spare1
       [-1   ,1   ], # 11Spare2
       [0    ,5   ], # 12Std of Speed
       [0    ,0.2 ], # 13Std of Intensity
       [0    ,1   ], # 14Track Density
       [-1   ,1   ], # 15X-component of Mean Orientation Vector
       [-40  ,5   ], # 16X-component of Mean Velocity
       [-1   ,1   ], # 17Y-component of Mean Orientation Vector
       [-40  ,5  ]] # 18Y-component of Mean Velocity

draw_var = ["fden","gden","lden","marea","mgdr","",
            "mlif","msp" ,"mstr","mten" ,""    ,"",
            ""    ,""    ,"tden",""    ] # 7 variables
#draw=[6,7,8]
draw=[14,1,2]
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
    flats = int(sys.argv[5])
    flatn = int(sys.argv[6])
    flonl = int(sys.argv[7])
    flonr = int(sys.argv[8])

if level == 1: #middle range 
    cnlvl[ 1][:]=[0    ,0.3 ]
    cnlvl[ 2][:]=[0    ,0.3 ]
    cnlvl[14][:]=[0    ,1.5 ]
if level == 2: # large range,use total cnlvlel bar 
    cnlvl[ 0][:]=[0   ,100 ]
    cnlvl[ 1][:]=[0.5 ,0.5 ]
    cnlvl[ 2][:]=[0.5 ,0.5 ]
    cnlvl[14][:]=[1   ,2   ]

figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig"
path = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev/statistic'
files = '%s/ff_250_1980-2020_stat.nc'%path
f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
del f,lat,lon

ds = xr.open_dataset('/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc')
da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
# increased performance by loading data into memory first, e.g., with load()
uwnd = da.mean('time')
del ds, da

ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

fig = plt.figure(figsize=(12,12),dpi=150)
ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
for nl in range(0,len(lev),1):
    files = '%s/%s_%d_1980-2020%s_stat.nc'%(path,prefix,lev[nl],suffix)
    print(files)
    f = xr.open_dataset(files)
    for nv in range(0,len(draw),1):#,len(f),1):
        var1 = f[draw_var[draw[nv]]][0,1,2]
        var = f[draw_var[draw[nv]]].sel(long=ilon,lat=ilat).mean("time").data
        print(var)
        if draw[nv] == 9:
            var=var*24
        if draw[nv] > 2 and draw[nv] != 14:
            tden = f['tden'].sel(long=ilon,lat=ilat).mean("time").data
            mask = tden < 1.0
            var = np.ma.array(var,mask=mask)
        if lev[nl]==850:
            var = np.ma.array(var,mask=(phis>1500))
        
        if cnlvl[draw[nv]][0] < 0 :
            fcolors = cmaps.BlueDarkRed18
        else:
            colr = cmaps.topo_15lev(range(0,16,1))[::-1]
            colr[8,:] = colors.to_rgba('y')
            colr[7,:] = colors.to_rgba('c')
            fcolors = colors.ListedColormap(colr)
            #fcolors = cmaps.precip2_17lev
        cnlevels = np.arange(cnlvl[draw[nv]][0], cnlvl[draw[nv]][0]+cnlvl[draw[nv]][1]*(
            fcolors.N-1), cnlvl[draw[nv]][1])
        norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
        
        axe = ax[nl][nv]
        axe.add_feature(cfeat.GSHHSFeature(scale='coarse',levels=[1,2],
            edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title("(%s) %d %s"%(numod[nl*3+nv],lev[nl],var1.long_name),fontdict=font)

        shad = axe.contourf(ilon, ilat, var, cnlevels, 
                     transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [500,1000,1500,3000,4500], 
                     transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        if dbox >= 1 :
            axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                     linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box

        jets = axe.contour(ilon, ilat, uwnd, [25,35,45,100], 
                     transform=ccrs.PlateCarree(),colors='red',linewidths=2)

        if nv == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nl == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
            position = fig.add_axes([xbar[nv], bmlo+0.05, 0.28, 0.008]) #left, bottom, width, height
            cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)

plt.tight_layout(rect=(0,bmlo,1,1))
plt.savefig("%s/stat_annual%s.png"%(figdir,suffix), bbox_inches='tight',pad_inches=0.01)


