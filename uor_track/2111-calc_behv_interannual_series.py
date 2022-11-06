#!/usr/bin/env python
'''
1. when cyclone reach the defined region, read datetime
    , one cyclone only one time be read
2. calc the cyclone number time series ([len(lev),len(behv),len(months),len(years)] 
3. plot the time series and store the data

2021-11-08
renql
'''

import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import pandas as pd
import sys, os, subprocess, linecache, gc
from datetime import datetime
from scipy import stats
from renql import cyc_filter
import matplotlib.pyplot as plt
from matplotlib import colors

def standardize(data):
    naxis = 0
    mu = np.mean(data, axis=naxis)
    sigma = np.std(data, axis=naxis)
    if mu == 0 and sigma == 0:
        data1 = np.zeros(data.shape)
    else:
        data1 = (data - mu) / sigma
    return data1 

def read_cyclone_time(filname,flats,flatn,flonl,flonr):
    ff = open(filname,"r") 
    line1 = ff.readline()
    line2 = ff.readline()
    line3 = ff.readline()
    line4 = ff.readline()
    a = line4.strip().split(" ",1)
    term2 = a[1].strip().split(" ",1)
    
    ct1=[]
    line = ff.readline()
    while line:
        term = line.strip().split(" ")
        if term[0] == "TRACK_ID":
            linenum = ff.readline()
            term1 =linenum.strip().split(" ")
            num = int(term1[-1])
            
            data=[]
            for nl in range(0,num,1):
                line = ff.readline()
                data.append(list(map(float,line.strip().split(" "))))
            
            for nl in range(0,num-1,1):
                if data[nl][1] <= flonl and data[nl+1][1] > flonl and data[nl+1][1] <= 180 :
                    point = cyc_filter.intersection_point_fixx([data[nl][1:3], data[nl+1][1:3]], flonl)
                    if point <= flatn and point >= flats :
                        a = str(int(data[nl+1][0]))
                        ct1.append(datetime.strptime(a,'%Y%m%d%H'))
                        break 
        line = ff.readline()
    ff.close()
    print("total cyclone number in %s : %s ; time number = %d" %(ff.name, term2[0], len(ct1)))
    return ct1

def calcbehv():
    var  = np.zeros( [len(lev),len(behv),len(months),len(year)],dtype=float )
    for nl in range(0,len(lev),1):
        for nr in range(0,len(behv),1):
            if nr == 0:
                filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix
            else:
                filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix+"_"+behv[nr]
            
            ct = read_cyclone_time(filname,flats,flatn,flonl,flonr)
            if len(ct) == 0:
                continue

            ctda = xr.DataArray(ct, coords=[ct], dims=["time"])
            for nm in range(0,len(months),1):
                if nm == 0 :
                    for ny in range(0,len(year),1):
                        index1 = np.array([ctda.time.dt.month.isin(12),ctda.time.dt.year.isin(year[ny]-1)]).all(axis=0)
                        index2 = np.array([ctda.time.dt.month.isin([1,2]),ctda.time.dt.year.isin(year[ny])]).all(axis=0)
                        ctda0  = ctda.sel(time=np.array([index1,index2]).any(axis=0))
                        var[nl,nr,nm,ny] = len(ctda0)
                else:
                    ctda0 = ctda.sel(time=ctda.time.dt.month.isin(months[nm]))
                    if len(ctda0) == 0:
                        continue
                    for ny in range(0,len(year),1):
                        ctda1 = ctda0.sel(time=ctda0.time.dt.year.isin(year[ny]))
                        var[nl,nr,nm,ny] = len(ctda1)
    return var

def regCoef_n(ts1,ts2,dim1,dim2):
    dims1 = ts1.shape
    dims2 = ts2.shape
    slope = np.zeros([dims1[0],dims1[1],dims2[1],dims2[2],dims2[3]],dtype=float)
    pvalu = np.zeros([dims1[0],dims1[1],dims2[1],dims2[2],dims2[3]],dtype=float)
   
    dims3 = slope.shape
    for nx0 in range(0,dims3[0]):
        for nx1 in range(0,dims3[1]):
            for nx2 in range(0,dims3[2]):
                for nx3 in range(0,dims3[3]):
                    for nx4 in range(0,dims3[4]):
                        regression = stats.linregress(ts1[nx0,nx1,:],ts2[:,nx2,nx3,nx4])
                        slope[nx0,nx1,nx2,nx3,nx4] = regression.slope
                        pvalu[nx0,nx1,nx2,nx3,nx4] = regression.pvalue
                        #print(regression.slope,regression.pvalue)
    return slope, pvalu

def new_linregress(x, y):
    # Wrapper around scipy linregress to use in apply_ufunc
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    return np.array([slope, p_value])

def calc_regression(var,varname,fileout):
    levc = [850,500,250,200]
    lonl=0  #0  #
    lonr=150#360#
    lats=0  #
    latn=90 #
    datapath = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/"#t/ERA5_NH_t_1989.nc
    ds  = xr.open_dataset(datapath+"t/ERA5_NH_t_1989.nc")
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=lonl) & (lon<=lonr)]
    ilat = lat[(lat>=lats) & (lat<=latn)]
    
    regr = np.zeros( [var.shape[0],var.shape[1],var.shape[2],len(levc),len(ilat),len(ilon)],dtype=float ) 
    prob = np.zeros( [var.shape[0],var.shape[1],var.shape[2],len(levc),len(ilat),len(ilon)],dtype=float ) 
    for nm in range(0,len(months),1):
        term = np.zeros( [len(year),len(levc),len(ilat),len(ilon)],dtype=float ) 
        for ny in range(0,len(year),1):
            print("ERA5_NH_%s_%d.nc"%(varname,year[ny]))
            ds   = xr.open_dataset("%s%s/ERA5_NH_%s_%d.nc"%(datapath,varname,varname,year[ny]))
            term[ny,:,:,:] = ds[varname].sel(time=ds.time.dt.month.isin(months[nm]),
                    level=levc,longitude=ilon,latitude=ilat).mean("time")

        foo = xr.DataArray(term, coords=[("year",year), ("lev",levc), ("lat",ilat), ("lon",ilon)]) 
        a1 = xr.apply_ufunc(new_linregress, var[:,:,nm,:],foo
                            input_core_dims=[['year'], ['year']],
                            output_core_dims=[["parameter"]],
                            vectorize=True,
                            dask="parallelized",
                            output_dtypes=['float64'],
                           )
        regr[:,:,nm,:,:,:] = a1[0,:,:,:,:,:] 
        prob[:,:,nm,:,:,:] = a1[1,:,:,:,:,:]

    ds = xr.Dataset(
            {
                "regr": (["lev3", "behv","month", "level", "lat","lon"], regr),
                "prob": (["lev3", "behv","month", "level", "lat","lon"], prob),
                },
            coords={
                "lev3" : (["lev3"],lev),
                "behv" : (["behv"],behv),
                "month": (["month"],season),
                "level": (["level"],levc),
                "lat"  : (["lat"],ilat.data),
                "lon"  : (["lon"],ilon.data),
                },
            )
    ds.attrs["description"] = "interannual regression "+suffix
    ds.to_netcdf(fileout+varname+".nc","w")
    return regr, prob

def storedata(var,fileout):
    ds = xr.Dataset(
            {
                "var": (["lev3", "behv","month","year"], var),
                },
            coords={
                "lev3" : (["lev3"],lev),
                "behv" : (["behv"],behv),
                "month": (["month"],season),
                "year" : (["year"],year),
                },
            )
    ds.attrs["description"] = '''annual time series for diff behavior of 3levs cyclones come from
        %d-%dE, %d-%dN, threshold timelife %d'''%(flonl,flonr,flats,flatn,time)
    ds.to_netcdf(fileout,"w")
    
def draw_ts(var,figdir):
# draw the interannual time series
    nrow = 4
    ncol = 1
    bmlo = 0.1
    bigfont=22
    midfont=14
    smfont=10
    plinsty = ["solid",(0,(3,2)),(0,(1,2))] # change with lev
    pcolor  = ["k","g","r","b"] # change with month 

    for nm in range(0,len(months),1):
        fig = plt.figure(figsize=(9,9),dpi=200)
        ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
        for npp in range(0,nrow,1):
            nb = npp+1
            ax[npp].set_title(season[nm]+" "+behv[nb]+" "+suffix,fontsize=midfont)
            for nl in range(0,len(lev),1):
                ax[npp].plot(year,var[nl,nb,nm,:],linewidth=1.2,color=pcolor[nl],
                        label=str(lev[nl])+": "+str(np.std(var[nl,nb,nm,:]))) #linestyle=plinsty[nl],
            ax[npp].set_xlim(year[0],year[-1])
       
        handles, labels = ax[0].get_legend_handles_labels()
        fig.legend(handles, labels, bbox_to_anchor=(0.1, bmlo-0.1, 0.8, 0.1),ncol=3)
        fig.tight_layout(rect=(0,bmlo,1,1)) #w_pad=0.5,h_pad=0.001) #,
        #fig.savefig(figdir)
        fig.savefig(figdir+season[nm]+".png", bbox_inches='tight',pad_inches=0.01)

def draw_relation(var,figdir):
    nrow = 2
    ncol = 2
    bmlo = 0.1
    bigfont=22
    midfont=10
    smfont=6
    cnlevels = [-0.99,-0.95,-0.90,0.90,0.95,0.99]
    listcolor = ["dodgerblue","deepskyblue","powderblue","white","gold","darkorange","red"]
    fcolors = colors.ListedColormap(listcolor)
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
    
    dims = var.shape
    var1 = var.reshape((dims[0]*dims[1],dims[2],dims[3]),order='C')
    
    label=[]
    for i0 in range(dims[0]):
        for i1 in range(dims[1]):
            label.append("%d %s"%(lev[i0],behv[i1+1]))

    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for ix in range(0,nrow,1):
        for iy in range(0,ncol,1):
            nm = 2*iy+ix
            cor  = np.empty((dims[0]*dims[1],dims[0]*dims[1]),dtype=float) 
            prob = np.ones((dims[0]*dims[1],dims[0]*dims[1]),dtype=float)
            for i in range(len(label)):
                for j in range(i):
                    cor[i,j], prob[i,j] = stats.pearsonr(var1[i,nm,:],var1[j,nm,:])
                    print("%s %s %s : r = %.2f prob = %.2f"%(season[nm], label[i], label[j], cor[i,j], prob[i,j]))
            prob = 1-prob
            prob = np.where(cor<0,-prob,prob)

            axe = ax[ix][iy]
            axe.set_title(season[nm]+" "+suffix,fontsize=midfont)
            im = axe.imshow(prob, cmap=fcolors, norm=norm)
            axe.set_xticks(np.arange(len(label)))
            axe.set_yticks(np.arange(len(label)))
            axe.set_xticklabels(label,fontsize=smfont)
            axe.set_yticklabels(label,fontsize=smfont)
            axe.plot([3.5,3.5],[0,(len(label)-1)],color="k")
            axe.plot([7.5,7.5],[0,(len(label)-1)],color="k")
            axe.plot([0,(len(label)-1)],[7.5,7.5],color="k")
            axe.plot([0,(len(label)-1)],[3.5,3.5],color="k")
            # Rotate the tick labels and set their alignment.
            plt.setp(axe.get_xticklabels(), rotation=45, ha="right",
                             rotation_mode="anchor")

            #axe.spines[:].set_visible(False)
            # Loop over data dimensions and create text annotations.
            for i in range(len(label)):
                for j in range(i):
                    text = axe.text(j,i,"%.2f"%cor[i,j], ha="center", va="center", color="k",fontsize=6)

    position = fig.add_axes([0.92, 0.2, 0.01, 0.6]) #left, bottom, width, height
    cb = plt.colorbar(im, cax=position ,orientation='vertical')#, shrink=.9)
    #fig.tight_layout(rect=(0,bmlo,1,1)) #w_pad=0.5,h_pad=0.001) #,
    #fig.savefig(figdir)
    fig.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.08)
            
if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 27 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 90 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "fft"
    season = 0 # 0 monthly, 1 seasonal
else:
    option= int(sys.argv[1]) 
    flats = int(sys.argv[2])
    flatn = int(sys.argv[3])
    flonl = int(sys.argv[4])
    flonr = int(sys.argv[5])
    prefix = int(sys.argv[6])
    season = int(sys.argv[7])
    time = int(sys.argv[8])

suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonl)
outdir="/home/users/qd201969/uor_track/mdata/"
figdir ="/home/users/qd201969/uor_track/fig/behv_ts_"
path = '/home/users/qd201969/ERA5-1HR-lev/'

behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
lev  = [850,500,250]
year = range(1980,2021,1)
months = [[12,1,2],[3,4,5],[6,7,8],[9,10,11]]
season = ["DJF","MAM","JJA","SON"]

#var = calcbehv()
#storedata(var,outdir+"behv_interannual_series.nc")

fvar = xr.open_dataset(outdir+"behv_interannual_series.nc")
var = fvar['var']

#draw_ts(var,figdir)
draw_relation(var[:,1:,:,:],figdir+"relation")

#for varname in ["u","v","z"]:
#    calc_regression(var,varname,outdir+"regr_interannual_")
#draw_regression()

