#!/usr/bin/env python
# monthly download hourly surface variable

import cdsapi
import numpy as np
import subprocess, os

c = cdsapi.Client()

def down_single_level(var,year):
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'grib',
            'expver': '1',
            'variable': [
                '10m_u_component_of_wind',
                '10m_v_component_of_wind',
            ],
            'year': str(year),
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': [
                '00:00', '01:00', '02:00',
                '03:00', '04:00', '05:00',
                '06:00', '07:00', '08:00',
                '09:00', '10:00', '11:00',
                '12:00', '13:00', '14:00',
                '15:00', '16:00', '17:00',
                '18:00', '19:00', '20:00',
                '21:00', '22:00', '23:00',
            ],
            'grid': [0.25, 0.25],
            'area': [70, 0, 15, 150, ],
        },
        '%s%s/ERA5_%s_%d.grib'%(path,var,var,year))

variables = {'precip':'total_precipitation','slp':'mean sea level pressure',
        'wind10':'10m u-component of wind','v10':'10m v-component of wind'}
path = '/work/scratch-pw2/renql/ERA5_hourly/'
#path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/'

for var in ('wind10',):
    for year in range(1979,2022):
        if os.path.isfile('%s%s/ERA5_%s_%d.grib'%(path,var,var,year)):
            print('%s%s/ERA5_%s_%d.grib exist'%(path,var,var,year))
        else:
            down_single_level(var,year)
#        com = "cdo expr,’speed=sqrt(sqr(uwnd)+sqr(vwnd))’ infile outfile"
#
#    os.chdir("%s%s"%(path,var))
#    for year in range(1981,2021):
#        infile1 = "ERA5_speed_%s_%d-12.nc"%(var,year-1)
#        infile2 = "ERA5_speed_%s_%d.nc"%(var,year)
#        infile3 = "ERA5_speed_%s_%d-1.nc"%(var,year+1)
#        com = "cdo -b F32 -mergetime %s %s %s ERA5_%s_1hr_dec-jan%d.nc"\
#            %(infile1,infile2,infile3,var,year)
#        ret=subprocess.Popen(com,shell=True)
#        ret.wait()
        #com = "rm ERA5_%s_%d-*.nc"%(var,year-1)
        #ret=subprocess.Popen(com,shell=True)
        #ret.wait()

