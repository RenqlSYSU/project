#!/usr/bin/env python
import cdsapi

varname = ['u_component_of_wind','v_component_of_wind','vorticity',
    'potential_vorticity','vertical_velocity','specific_humidity',]
filname = ['u','v','vo','pv','omega','q']
path = '/gws/nopw/j04/ncas_generic/users/renql/ERA5_mon'

c = cdsapi.Client()

def main_run():
#    for nv in range(3,len(varname),1):
#        down_pressure_level(varname[nv],filname[nv])
    #down_single_level('2m_temperature','t2')    
    down_pressure_level('vorticity','vo')

def down_single_level(var,filvar):
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': var,
            'year': [
                '1979', '1980', '1981',
                '1982', '1983', '1984',
                '1985', '1986', '1987',
                '1988', '1989', '1990',
                '1991', '1992', '1993',
                '1994', '1995', '1996',
                '1997', '1998', '1999',
                '2000', '2001', '2002',
                '2003', '2004', '2005',
                '2006', '2007', '2008',
                '2009', '2010', '2011',
                '2012', '2013', '2014',
                '2015', '2016', '2017',
                '2018', '2019', '2020',
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
            'expver': '1',
        },
        '%s/ERA5_mon_%s_1979-2020.nc'%(path,filvar))

def down_pressure_level(var,filvar):
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'format': 'netcdf',
            'product_type': 'monthly_averaged_reanalysis',
            'variable': var,
            'year': [
                '1979', '1980', '1981',
                '1982', '1983', '1984',
                '1985', '1986', '1987',
                '1988', '1989', '1990',
                '1991', '1992', '1993',
                '1994', '1995', '1996',
                '1997', '1998', '1999',
                '2000', '2001', '2002',
                '2003', '2004', '2005',
                '2006', '2007', '2008',
                '2009', '2010', '2011',
                '2012', '2013', '2014',
                '2015', '2016', '2017',
                '2018', '2019', '2020',
            ],
            'month': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
            ],
            'time': '00:00',
            'expver': '1',
            'pressure_level': [
                '200', '250', '300',
                '500', '850',
            ],
        },
        '%s/ERA5_mon_%s_1979-2020.nc'%(path,filvar))

if __name__=='__main__':
    main_run()

