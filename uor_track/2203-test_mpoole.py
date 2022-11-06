#!/usr/bin/env python
import xarray as xr
import numpy as np
import pandas as pd
from multiprocessing import Pool,Manager

def test(i,temp):
    var=temp.get()
    print('task %d, var[0] before change %f'%(i,var[0,i]))
    #print(i,var)
    for j in range(len(var[:,i])):
        var[j,i] = var[j,i] + 1
    print('task %d, var[0] after change %f'%(i,var[0,i]))
    temp.put(var)
    #return var

data = np.random.rand(4, 3)
locs = ["IA", "IL", "IN"]
times = pd.date_range("2000-01-01", periods=4)
foo = xr.DataArray(data, coords=[times, locs], dims=["time", "space"])
print(foo)

temp = Manager().Queue()
temp.put(foo)
print(temp)

process_pool = Pool(processes=3)
results=[]
for nl in range(3):
    result = process_pool.apply_async(test, args=(nl,temp))
    results.append(result)
print(results) 
print(results[0].get()) 
process_pool.close()
process_pool.join() 
print("finish")
#for nl in range(3):
#    foo[:,nl] = results[nl].get()
print(foo)
print(temp)
print(temp.get())

