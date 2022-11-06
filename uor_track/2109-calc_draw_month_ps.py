#!/usr/bin/env python
import xarray as xr
import numpy as np
import subprocess
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30 #60 #

path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon/'
figdir = '/home/users/qd201969/uor_track/fig'
var_name = ['u','sp']
nv = 1

# -------------- read data and calc climatology ---------------
ds = xr.open_dataset(path+'ERA5_mon_'+var_name[nv]+'_1979-2020.nc')
lat = ds.latitude
lon = ds.longitude
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
if nv == 1:
    da = ds[var_name[nv]].sel(longitude=ilon,latitude=ilat)
    da /= 100.0 # convert Pa to hPa
    lev  =[400,1200,50]
else:
    da = ds[var_name[nv]].sel(level=200,expver=1,longitude=ilon,latitude=ilat)
    lev  =[20,52,2]
# increased performance by loading data into memory first, e.g., with load()
var = da.groupby(da.time.dt.month).mean('time').data

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").data
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

# -------------- draw figure ---------------
titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
cnlevels = np.array([450, 500, 550, 825, 850, 875, 900, 950])
fcolors = colors.ListedColormap(['white','k','c','m','y', "coral","blue","green",'pink'])
norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

nrow = 4 #6 #
ncol = 3 #2 #
bmlo = 0.35 #0.25 #
fig = plt.figure(figsize=(12,12),dpi=300)
ax = fig.subplots(nrow, ncol, subplot_kw=dict(
    projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
nm = -1
for nr in range(0,nrow,1):
    for nc in range(0,ncol,1):
        nm = nm+1 
        axe = ax[nr][nc]
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title(titls[nm],fontsize=title_font)

        cont = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                     transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                     transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
        if nc == 0:
            axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
            axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nr == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
plt.savefig('%s/%s.png'%(figdir,var_name[nv]), bbox_inches='tight',pad_inches=0.01)

