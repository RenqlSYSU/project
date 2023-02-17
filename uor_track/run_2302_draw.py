#!/usr/bin/env python
import sys, os, subprocess

suffixs = ['_6local','_6outside']
path = "/home/ys17-23/Extension2/renql/project/uor_track/"

for suffix in suffixs: 
    #com = "python %s/2302-draw_weather_contri.py '%s' contri"%(path,suffix)
    com = "python %s/2302-draw_weather_contri.py '%s' maxpreci"%(path,suffix)
    print(com)
    ret=subprocess.Popen(com,shell=True)
    ret.wait()

