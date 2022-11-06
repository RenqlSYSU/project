#!/usr/bin/env python
'''
read total cyclone number in ff_250_1980-2020_2_3045-6060
then judge whether the cyclone is NTN, STN or LYS in this region
if not, it is considered to be able to pass over the TP
plot the filter box in the trajectory figure

20211014
renql
'''

import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import sys, os, subprocess, linecache, gc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from renql import cyc_filter, monthly_calc, life_intensity

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 27 #int(sys.argv[2])
    flatn = 45 #int(sys.argv[3])
    flonl = 60 #int(sys.argv[4])
    flonr = 85 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "ff"
    suffix0 = "_5_2545-60110"
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

suffix=str(option)+"_"+str(flatn)+str(flatn)+"-"+str(flonl)+str(flonr)
figdir = "/home/users/qd201969/uor_track/fig/"
fileout="/home/users/qd201969/uor_track/mdata/behv3_month_%dh_%s.nc"%(time,suffix)
calcbehv = 1
draw_hist = 0
storedata  = 0
drawannual = 0
drawbox    = 0

flonr2 = 90
behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
lats = [flatn ,flatn ,flats ,flats ,flats ]
latn = [flatn ,flatn ,flats ,flatn ,flatn ]
lonl = [flonl ,flonl ,flonl ,flonr2,flonl ]
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
path = '/home/users/qd201969/ERA5-1HR-lev/'
os.chdir("/home/users/qd201969/uor_track/fig")
#===============================================
# box filer cyclone, read cyclone number
#===================================================
if calcbehv == 1:
    var  = np.empty( [len(month),len(lev),len(behv)],dtype=int )  
    perc = np.empty( [len(month),len(lev),len(behv)],dtype=float )  
    for nl in range(0,len(lev),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0+"_"+suffix#+"_"+str(time)

        if not os.path.isfile(filname) :
            filname0 = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0
            if option == 2:
                var[0,nl,0] = cyc_filter.line_filt(filname0,flatn,flatn,flonl,flonr,option,time)
            else:
                com = "/home/users/qd201969/track-1.5.2/utils/bin/box %s %d %d %d %d %d 0 0.0"\
                    %(filname0,flats,flatn,flonl,flonr,option)
                ret=subprocess.popen(com,shell=True)
                ret.wait()
                subprocess.run("mv %s.new %s"%(filname0,filname),shell=True)

        if not os.path.isfile(filname+"_NTN") :
            var[0,nl,:] = cyc_filter.behavior(filname,
                    lats[1:len(behv)],latn[1:len(behv)],lonl[1:len(behv)],lonr[1:len(behv)])

        for nr in range(0,len(behv),1):
            if nr == 0:
                filname1 = filname
            else:
                filname1 = filname+"_"+behv[nr]

            #life_intensity.hist_life_intensity(filname1)

            #var[1:len(month),nl,nr] = monthly_calc.calc_month_fft(filname1,time)
            #var[1:len(month),nl,nr] = monthly_calc.calc_month(
            #        frae, month[1:len(month)], nday[1:len(month)], filname1)
            var[1:len(month),nl,nr] = monthly_calc.draw_month_traj(
               frae, month[1:len(month)],nday[1:len(month)],filname1,lats,latn,lonl,lonr)

if draw_hist == 1:
    nrow = len(lev)
    ncol = len(behv)-1
    bmlo = 0.1
    BIGFONT=22
    MIDFONT=14
    SMFONT=10
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol, sharex=True, sharey=True)
    for nl in range(0,len(lev),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0+"_"+suffix#+"_"+str(time)
        for nr in range(1,len(behv),1):
            if nr == 0:
                filname1 = filname
            else:
                filname1 = filname+"_"+behv[nr]
            
            patches = life_intensity.hist_life_intensity(filname1, \
                    ax=ax[nl,nr-1], title="%d %s"%(lev[nl],behv[nr]))
            if nl == 2:
                ax[nl,nr-1].set_xlabel("lifetime (days)",fontsize=SMFONT)

            if nr == 1:
                ax[nl,nr-1].set_ylabel("Max intensity ($10^{-5} s^{-1}$)",fontsize=SMFONT)

    position = fig.add_axes([0.96, bmlo+0.05, 0.02, 0.8]) #left, bottom, width, height
    cb = plt.colorbar(patches, cax=position ,orientation='vertical')#, shrink=.9)
    fig.tight_layout(rect=(0,bmlo,0.94,1)) #w_pad=0.5,h_pad=0.001) #,
    fig.savefig("%slife_inte_%s.png"%(figdir,suffix), bbox_inches='tight',pad_inches=0.08)

if storedata == 1:
    f = open("%sbehv3_month_%dh_%s"%(figdir,time,suffix),'w')
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
    
    ax[2][1].legend(behv[1:len(behv)], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sbehv3_month_%dh_%s"%(figdir,time,suffix))
    #fig.savefig("%sbehv3_month_%dh_%s"%(figdir,time,suffix), bbox_inches='tight',pad_inches=0.01)
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax  = fig.subplots(1, 1) #sharex=True, sharey=True
    ax.set_title("total number "+suffix,fontsize=midfont)
    for nl in range(0,len(lev),1):
        ax.plot(imonth,var[1:len(month),nl,0],linewidth=3,label=str(lev[nl]))
    
    ax.legend(loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sbehv3_month_%dh_%s_total.png"%(figdir,time,suffix))

#===============================================
# draw trajectory and filter box 
#===================================================
if drawbox == 1:
    if calcbehv != 1:
        fvar = xr.open_dataset(fileout)
        var = fvar['num'].data
        perc= fvar['perc'].data

    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    cfp.setvars(file="traj-"+prefix+"_"+suffix+".png")
    cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.02,bottom=0.5) #0.55
    cfp.mapset(lonmin=0, lonmax=150, latmin=15, latmax=70)
    for nr in range(1,len(lats),1):
        if nr == 0:
            suffix2 = ""
        else:
            suffix2 = "_"+behv[nr]
        
        for nl in range(0,3,1):#,len(f),1):
            np=(nr-1)*3+nl+1
            
            cfp.gpos(np)
            cfp.levs(manual=[1500,3000,4500])
            cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,\
                    title=suffix+" "+str(lev[nl])+" "+behv[nr]+" "+
                    str(var[0,nl,nr])+" "+str(round(perc[0,nl,nr],1))+"%")
            cfp.levs()
            
            for nrr in range(0,len(lats),1):
                cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],
                        [latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                        linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box

            if var[0,nl,nr] == 0:
                continue
            
            filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0+"_"+suffix+suffix2
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
            
    cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title='Relative Vorticity (Hz)*1e5') #0.53
    cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./traj-"+prefix+"_"+suffix+".png",shell=True) 
    #subprocess.run("rm /home/users/qd201969/ERA5-1HR-lev/*.nc",shell=True) 

