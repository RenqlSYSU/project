#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-6060
then use 2110-line_filter.py to filter different behaviors cyclones
plot the filter box in the trajectory figure

20210928
'''

import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import sys, os, subprocess, linecache, gc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def calc_month( frae, month, nday, filname):
    var  = np.empty( [len(month)], dtype=int ) # 4 or 12 
    datafile = "/home/users/qd201969/uor_track/fig/match_ens1_yes.dat"
    for nm in range(0,len(month),1):
        fras = frae+1
        frae = fras+24*nday[nm]-1
        print("month=%s, frame_s=%d, frame_e=%d" %(month[nm],fras,frae))
        rf = open("/home/users/qd201969/uor_track/fig/region.dat","w")
        rf.write("0 360\n-90 90\n"+str(fras)+" "+str(frae))
        rf.close()

        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/censemble2 "+\
                filname+" "+filname+" 0 100 10 1 0 0 0 0 s 0 1 > rec",shell=True)
        ret.wait()

        a=linecache.getline(datafile, 4).strip().split(" ",1)
        term = a[1].strip().split(" ",1)
        var[nm] = int(term[0])
        linecache.clearcache()
    return var

def draw_month_traj(frae, month, nday, filname, behv):
    var  = np.empty( [len(month)], dtype=int ) # 4 or 12
    fterm = filname.split("_")
    datafile = "/home/users/qd201969/uor_track/fig/match_ens1_yes.dat"
    figname = "traj-mon_"+fterm[1]+"_"+fterm[4]+"_"+behv+".png"
    
    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    phis=phis/9.8 # transfer from m2/s2 to m

    mlonl=0  #0  #
    mlonr=150#360#
    mlats=15 #0  #
    mlatn=70 #90 #
    ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
    lat = ds.latitude
    lon = ds.longitude
    ilon = lon[(lon>=mlonl) & (lon<=mlonr)]
    ilat = lat[(lat>=mlats) & (lat<=mlatn)]
    da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat).load()
    # increased performance by loading data into memory first, e.g., with load()
    uwnd = da.groupby(da.time.dt.month).mean('time')
    del ds, da
    gc.collect()
   
    cfp.setvars(file=figname)
    cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
    cfp.mapset(lonmin=mlonl, lonmax=mlonr, latmin=mlats, latmax=mlatn)
    for nm in range(0,len(month),1):
        fras = frae+1
        frae = fras+24*nday[nm]-1
        rf = open("/home/users/qd201969/uor_track/fig/region.dat","w")
        rf.write("0 360\n-90 90\n"+str(fras)+" "+str(frae))
        rf.close()
        
        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/censemble2 "+\
                filname+" "+filname+" 0 100 10 1 0 0 0 0 s 0 1 > rec",shell=True)
        ret.wait()
        a=linecache.getline(datafile, 4).strip().split(" ",1)
        term = a[1].strip().split(" ",1)
        var[nm] = int(term[0])
        linecache.clearcache()
        print("month=%s, frame_s=%d, frame_e=%d, number= %d" %(month[nm],fras,frae,var[nm]))
        
        cfp.gpos(nm+1)
        cfp.levs(manual=[1500,3000,4500])
        cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,line_labels=False,
                title=behv+" "+fterm[1]+"hPa "+fterm[4]+" "+month[nm]+" "+str(var[nm]))
        cfp.levs(manual=[30,32]) # jet stream
        cfp.con(uwnd[nm,:,:],x=ilon, y=ilat, ptype=1, fill=False,line_labels=False,
                lines=True, colors='indigo',linewidths=3.5)
        for nrr in range(0,len(lats),1):
            cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],
                    [latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                    linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box
        
        if var[nm] == 0:
            continue
        
        ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/tr2nc \
                "+datafile+" s /home/users/qd201969/TRACK-1.5.2/utils/TR2NC/tr2nc.meta",shell=True)
        ret.wait()
        f=cf.read(datafile+'.nc')
        g = f[2]
        g = g*1e5
        cfp.cscale('precip2_17lev')
        cfp.levs(min=0.0, max=8.0, step=0.5)
        cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
        
    cfp.cscale('precip2_17lev')
    cfp.levs(min=0.0, max=8.0, step=0.5)
    cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title='Relative Vorticity (Hz)*1e5')
    cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./"+figname,shell=True) 
    return var

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 25 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 60 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "ff"
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
suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
figdir = "/home/users/qd201969/uor_track/fig/behv_month_"+str(time)+"h_"+suffix
fileout="/home/users/qd201969/uor_track/mdata/behv_month_"+str(time)+"h_"+suffix+".nc"
calcbehv = 1
drawannual = 1
drawbox = 1

flonr2 = 90
behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
lats = [flats ,flatn ,flats ,flats ,flats ]
latn = [flatn ,flatn ,flats ,flatn ,flatn ]
lonl = [flonl ,flonl ,flonl ,flonr2,flonr ]
lonr = [flonr ,flonr2,flonr2,flonr2,flonr2]
opti = [option,2     ,2     ,2     ,1     ]
lev = [850,500,250]

if season == 0:
    nday=[365,31,28,31,30,31,30,31,31,30,31,30,31]
    month=["ALL","Jan","Feb","Mar","Apr","May","Jun",\
            "Jul","Aug","Sep","Oct","Nov","Dec"]
    frae=744
else:
    nday=[365,90,92,92,91]
    month=["ALL","DJF","MAM","JJA","SON"]
    frae=0

imonth=range(1,len(month),1)
var  = np.empty( [len(month),len(lev),len(behv)],dtype=int )  
perc = np.empty( [len(month),len(lev),len(behv)],dtype=float )  
path = '/home/users/qd201969/ERA5-1HR-lev/'
os.chdir("/home/users/qd201969/uor_track/fig")
#===============================================
# box filer cyclone, read cyclone number
#===================================================
if calcbehv == 1:
    for nl in range(0,len(lev),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix#+"_"+str(time)

        if not os.path.isfile(filname) :
            filname0 = path+prefix+"_"+str(lev[nl])+"_1980-2020"
            ret=subprocess.Popen("python /home/users/qd201969/uor_track/2110-line_filter.py "+ \
                    filname0+" "+str(lats[0])+" "+str(latn[0])+\
                    " "+str(lonl[0])+" "+str(lonr[0])+" "+str(opti[0])+" "+str(time),shell=True)
            ret.wait()

        a=linecache.getline(filname, 4).strip() #读取特定行数据
        a=a.split(" ",1)
        term = a[1].strip().split(" ",1)
        var[0,nl,0] = int(term[0])
        linecache.clearcache()
        var[1:len(month),nl,0] = calc_month( frae, month[1:len(month)], nday[1:len(month)], filname)
        #var[1:len(month),nl,0] = draw_month_traj( frae, month[1:len(month)], nday[1:len(month)], filname, behv[0])

        for nr in range(1,len(lats)-1,1):
            suffix2 = str(opti[nr])+"_"+str(lats[nr])+str(latn[nr])+"-"+str(lonl[nr])+str(lonr[nr])
            if behv[nr] == "lys" :
                ret=subprocess.Popen("/home/users/qd201969/track-1.5.2/utils/bin/box "+\
                        filname+" "+str(lats[nr])+" "+str(latn[nr])+\
                        " "+str(lonl[nr])+" "+str(lonr[nr])+" "+str(opti[nr])+" 0 0.0",shell=True)
                ret.wait()
                subprocess.run("mv "+filname+".new "+filname+"_"+suffix2,shell=True)
            else:
                ret=subprocess.Popen("python /home/users/qd201969/uor_track/2110-line_filter.py "+ \
                        filname+" "+str(lats[nr])+" "+str(latn[nr])+\
                        " "+str(lonl[nr])+" "+str(lonr[nr])+" "+str(opti[nr])+" "+str(time),shell=True)
                ret.wait()
            
            a=linecache.getline(filname+"_"+suffix2, 4).strip().split(" ",1)
            term = a[1].strip().split(" ",1)
            var[0,nl,nr] = int(term[0])
            linecache.clearcache()
            var[1:len(month),nl,nr] = calc_month( frae, month[1:len(month)], nday[1:len(month)], filname+"_"+suffix2)
            #var[1:len(month),nl,nr] = draw_month_traj( frae, month[1:len(month)], nday[1:len(month)], filname+"_"+suffix2, behv[nr])

    var[:,:,-1] = var[:,:,0]-np.sum(var[:,:,1:4],axis=2)

    f = open("/home/users/qd201969/uor_track/behv-monthly-"+str(time)+"h_"+prefix+"-"+suffix,'w')
    f.write(path+prefix+"_1980-2020_"+suffix+"\n")
    f.write("*"*50+"\n")
    f.write("lats"+str(lats)+"\n")
    f.write("latn"+str(latn)+"\n")
    f.write("lonl"+str(lonl)+"\n")
    f.write("lonr"+str(lonr)+"\n")
    f.write("opti"+str(opti)+"\n")
    f.write("*"*50+"\n")
    for nm in range(0,len(month),1):
        f.write("\n")
        f.write("* "+month[nm]+" number (percent)\n")
        f.write("* lev, "+str(behv).strip('[').strip(']').replace('\'','')+"\n")
        for nl in range(0,len(lev),1):
            perc[nm,nl,:] = 100*var[nm,nl,:]/var[nm,nl,0]
            output = str(var[nm,nl,0])
            for nr in range(1,len(lats),1):
                output = output+", "+str(var[nm,nl,nr])+" ("+str(np.around(perc[nm,nl,nr],decimals=2))+")"
            f.write("* "+str(lev[nl])+", "+output+"\n")
    f.close()

    ds = xr.Dataset(
            {
                "num" : (["month", "lev", "behv"], var),
                "perc": (["month", "lev", "behv"], perc),
                },
            coords={
                "month": range(0,len(month),1), 
                "lev"  : (["lev"], lev),
                "behv" : (["behv"],behv),
                },
            )
    ds.attrs["description"] = "number and percent of cycle, minimum lifetime %d hour"%time
    ds.to_netcdf(fileout,"w")

#===============================================
# draw figure 
#===================================================
if drawannual == 1:
    nrow = 3
    ncol = 2
    bmlo = 0.4
    bigfont=22
    midfont=14
    smfont=10
    if calcbehv != 1:
        fvar = xr.open_dataset(fileout)
        var = fvar['num'].data
        perc= fvar['perc'].data

    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for nl in range(0,len(lev),1):
        ax[nl][0].set_title(str(lev[nl])+" number "+suffix,fontsize=midfont)
        ax[nl][1].set_title(str(lev[nl])+" percent "+suffix,fontsize=midfont)
        for nr in range(1,len(lats),1):
            ax[nl][0].plot(imonth,var[1:len(month),nl,nr],linewidth=3)
            ax[nl][1].plot(imonth,perc[1:len(month),nl,nr],linewidth=3)
    
    ax[2][1].legend(behv[1:5], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig(figdir+".png")
    #fig.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.01)

#===============================================
# draw trajectory and filter box 
#===================================================
if drawbox == 1:
    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    for nr in range(0,len(lats)-1,1):
        if nr == 0:
            suffix2 = ""
        else:
            suffix2 = "_"+str(opti[nr])+"_"+str(lats[nr])+str(latn[nr])+"-"+str(lonl[nr])+str(lonr[nr])
        
        cfp.setvars(file="traj-"+prefix+"_"+suffix+suffix2+".png")
        cfp.gopen(rows=3, columns=1 ,hspace=0.25)#,bottom=0.2
        cfp.mapset(lonmin=0, lonmax=150, latmin=10, latmax=70)
        for nl in range(0,3,1):#,len(f),1):
            np=nl+1
            
            cfp.gpos(np)
            cfp.levs(manual=[1500,3000,4500])
            cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,\
                    title=suffix+" "+str(lev[nl])+" "+behv[nr]+" "+
                    str(var[0,nl,nr])+" "+str(perc[0,nl,nr])+"%")
            cfp.levs()
            
            for nrr in range(0,len(lats),1):
                cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],
                        [latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                        linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box

            if behv[nr] == "PAS" and lev[nl] == 850:
                continue
            
            filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020_"+suffix+suffix2
            if not os.path.isfile(filname+'.nc') :
                ret=subprocess.Popen("/home/users/qd201969/TRACK-1.5.2/utils/bin/tr2nc \
                        "+filname+" s /home/users/qd201969/TRACK-1.5.2/utils/TR2NC/tr2nc.meta",shell=True)
                ret.wait()
            f=cf.read(filname+'.nc')
            print(f)
            g = f[2]
            g = g*1e5
            print(g)

            cfp.cscale('precip2_17lev')
            cfp.levs(min=0.0, max=8.0, step=0.5)
            cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
            
        cfp.cbar(position=[0.75, 0.2, 0.01, 0.6], title='Relative Vorticity (Hz)*1e5',orientation='vertical')
        cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./traj-"+prefix+"_"+suffix+"*.png",shell=True) 

