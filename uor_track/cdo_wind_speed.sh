#!/bin/bash

#cd /gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/z/
cd /work/scratch-pw2/renql/ERA5_hourly/wind10/ 

for ny in {1979..2021};do
    if [ -f ERA5_wind10_${ny}.grib ];then
        cdo expr,'speed=sqrt(sqr(var165)+sqr(var166))' \
            ERA5_wind10_${ny}.grib ERA5_speed10_${ny}.grib
        rm ERA5_wind10_${ny}.grib
        cdo -f nc copy ERA5_speed10_${ny}.grib ERA5_speed10_${ny}.nc
        rm ERA5_speed10_${ny}.grib 
    fi
done

for ny in {1980..2020};do
    ny1=$((ny+1))
    ny0=$((ny-1))
    if [ ! -f ERA5_speed10_1hr_dec-jan${ny}.nc -a -f ERA5_speed10_${ny1}.nc ];then
        cdo selmonth,1 ERA5_speed10_${ny1}.nc ERA5_speed10_${ny1}-01.nc
        cdo selmonth,12 ERA5_speed10_${ny0}.nc ERA5_speed10_${ny0}-12.nc
    
        cdo -b F32 -mergetime ERA5_speed10_${ny0}-12.nc ERA5_speed10_${ny}.nc \
            ERA5_speed10_${ny1}-01.nc ERA5_speed10_1hr_dec-jan${ny}.nc   
        
        rm ERA5_speed10_${ny1}-01.nc
        rm ERA5_speed10_${ny0}-12.nc
        rm ERA5_speed10_${ny0}.nc
    fi 
done

