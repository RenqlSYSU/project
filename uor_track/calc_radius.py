#!/usr/bin/env python
import numpy as np

def calc_radius(lat,wavenumber):
    a=6371 #km
    return a*np.cos(lat*np.pi/180.0)*2*np.pi/wavenumber/2

def calc_distance(lat,dlon):
    a=6371 #km
    return a*np.cos(lat*np.pi/180.0)*dlon*np.pi/180

print('20N, 7 : %f km'%calc_distance(15,7))
print('60N, 13 : %f km'%calc_distance(60,13))

'''
print('15N, 63 : %f km'%calc_radius(15,63))
print('70N, 63 : %f km'%calc_radius(70,63))
print('45N, 63 : %f km'%calc_radius(45,63))

print('15N, 5 : %f km'%calc_radius(15,5))
print('70N, 5 : %f km'%calc_radius(70,5))
'''
