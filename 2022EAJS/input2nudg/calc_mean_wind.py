import xarray as xr
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

path = '/WORK/sysu_hjkx_ys/renql/JET_NUDG/input2nudg'
time1 = datetime(1995,5,20).timetuple().tm_yday-1
time2 = datetime(1995,6,10).timetuple().tm_yday

ds = xr.open_dataset('%s/F2000_CAM5.cam.h1.ESM.clim.U.nc'%path)
da = ds['U'][time1:time2,:,:,:]
print(da.time)

var = da.mean('time')
print(var)

ds2 = var.to_dataset(name="U")
ds2.to_netcdf('%s/F2000_CAM5.cam.h1.ESM.clim.month_U.nc'%path,"w")

