import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
import cmaps

filein='/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.'
fileout='/home/ys17-19/renql/project/TP_NUDG/analysis/mdata/CTRL-Clim_daily_preci.nc'

lats = 15 
latn = 55 
lonl = 50  
lonr = 130

ds = xr.open_dataset('%s1979.PRECL.nc'%filein)
lat = ds.lat
lon = ds.lon
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
var=ds['PRECL'].sel(lon=ilon,lat=ilat)
ds = xr.open_dataset('%s1979.PRECC.nc'%filein)
var=ds['PRECL'].sel(lon=ilon,lat=ilat)

for ny in range(1979,2006):
print(var)
print(var.groupby(ds.time.dt.year))
#ds = xr.open_mfdataset('%s*PRECC.nc'%filein,concat_dim='time',combine='nested')
