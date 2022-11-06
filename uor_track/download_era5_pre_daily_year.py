#!/usr/bin/env python
# monthly download hourly surface variable

import cdsapi
import numpy as np
import subprocess, os

c = cdsapi.Client()

day28 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28',]

day29 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28','29',]

day30 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28', '29', '30',]

day31 = ['01', '02', '03','04', '05', '06', '07', '08', '09',
         '10', '11', '12', '13', '14', '15', '16', '17', '18',
         '19', '20', '21', '22', '23', '24', '25', '26', '27',
         '28', '29', '30', '31',]

mon_day    = {1:day31, 2:day28, 3:day31, 4:day30, 5:day31, 6:day30, 7:day31, 8:day31, 9:day30, 10:day31, 11:day30, 12:day31, }
mon_day_29 = {1:day31, 2:day29, 3:day31, 4:day30, 5:day31, 6:day30, 7:day31, 8:day31, 9:day30, 10:day31, 11:day30, 12:day31, }

year_mon = {}
for year in range(1979,2022):
    if ((year%4==0) and (year%100 !=0)) or (year%400)==0: 
        year_mon[year] = mon_day_29
    else:
        year_mon[year] = mon_day

def down_single_level(var,year,mon):
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'format': 'netcdf',
            'expver': '1',
            'variable': [
                '10m_u_component_of_wind',
                '10m_v_component_of_wind',
            ],
            'year': str(year),
            'month': [
                str(mon),
            ],
            'day': year_mon[year][mon],
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
        },
        '%s%s/ERA5_%s_%d-%d.nc'%(path,var,var,year,mon))

variables = {'precip':'total_precipitation','slp':'mean sea level pressure',
        'u10':'10m u-component of wind','v10':'10m v-component of wind'}
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_hourly/'

for var in ('precip',):
    #down_single_level(var,1979,12)
    #down_single_level(var,2021,1)
    #for year in range(1980,2021):
    #    for mon in range(1,13):
    #        down_single_level(var,year,mon)

    os.chdir("%s%s"%(path,var))
    for year in range(2001,2021):
        infile1 = "ERA5_%s_%d-12.nc"%(var,year-1)
        infile2 = "ERA5_%s_%d-[1-9].nc"%(var,year)
        infile3 = "ERA5_%s_%d-1[0-2].nc"%(var,year)
        infile4 = "ERA5_%s_%d-1.nc"%(var,year+1)
        com = "cdo -b F32 -mergetime %s %s %s %s ERA5_%s_1hr_dec-jan%d.nc"\
            %(infile1,infile2,infile3,infile4,var,year)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()
        com = "rm ERA5_%s_%d-*.nc"%(var,year-1)
        ret=subprocess.Popen(com,shell=True)
        ret.wait()

#'variable': variables[var],
