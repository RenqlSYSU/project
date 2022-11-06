#!/usr/bin/env python
'''
first line filter to get ff_250_1980-2020_2_3045-6060
then do the three level match using the shell script
then judge whether the cyclone is NTN, STN or PAS in this region
if not, it is considered to be lysis
plot the filter box in the trajectory figure

20211018
renql
'''

import cf
import cfplot as cfp
import xarray as xr
import numpy as np
import sys, os, subprocess, linecache, gc
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from renql import cyc_filter, monthly_calc

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 27 #int(sys.argv[2])
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
figdir = "/home/users/qd201969/uor_track/fig/behv3_match2_%dh_%s"%(time,suffix)
fileout="/home/users/qd201969/uor_track/mdata/behv3_match2_%dh_%s.nc"%(time,suffix)
calcbehv =0
drawannual = 0
drawbox = 1

flonr2 = 90
behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
lats = [flats ,flatn ,flats ,flats ,flats ]
latn = [flatn ,flatn ,flats ,flatn ,flatn ]
lonl = [flonl ,flonl ,flonl ,flonr2,flonr ]
lonr = [flonr ,flonr2,flonr2,flonr2,flonr2]
opti = [option,2     ,2     ,2     ,1     ]
lev1 = [850,500,250]
lev2 = ['only_250',          '250-500','250-500-850',
        'only_500','500-850','500-250','500-250-850', 
        'only_850','850-500',          '850-500-250']
npp = [10, 7, 1,  11, 5, 8, 2,  12, 6, 3]

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
var  = np.empty( [len(month),len(lev2),len(behv)],dtype=int )  
perc = np.empty( [len(month),len(lev2),len(behv)],dtype=float )  
perc2= np.empty( [len(month),len(lev2),len(behv)],dtype=float )  
path = '/home/users/qd201969/ERA5-1HR-lev/'
#===============================================
# box filer cyclone, read cyclone number
#===================================================
if calcbehv == 1:
    os.chdir("/home/users/qd201969/TRACK-1.5.2")
    for nl in range(0,len(lev1),1):
        filname  = path+prefix+"_"+str(lev1[nl])+"_1980-2020_"+suffix#+"_"+str(time)
        if not os.path.isfile(filname) :
            filname0 = path+prefix+"_"+str(lev1[nl])+"_1980-2020"
            var[0,nl,0] = cyc_filter.line_filt(filname0,flats,flatn,flonl,flonr,option,time)
    
    ret=subprocess.Popen("sh /home/users/qd201969/uor_track/match_3lev.sh %s _%s"
        %(prefix,suffix),shell=True)
    ret.wait()
    ff = open("%smatch_%s/number"%(path,suffix),'r')
    lines = ff.readlines()
    for nl in range(0,len(lines),1):
        term = lines[nl].strip().split(" ")
        var[0,nl,:] = cyc_filter.behavior(term[0].strip(),
                lats[1:len(behv)],latn[1:len(behv)],lonl[1:len(behv)],lonr[1:len(behv)])
        #var[1:len(month),nl,nr] = monthly_calc.calc_month(
        #        frae, month[1:len(month)], nday[1:len(month)],term[0].strip()) 
        #var[1:len(month),nl,nr] = monthly_calc.draw_month_traj(
        #   frae, month[1:len(month)],nday[1:len(month)],term[0].strip(),lats,latn,lonl,lonr)
    ff.close()
    del lines

    f = open(figdir,'w')
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
        f.write("* structure, "+str(behv).strip('[').strip(']').replace('\'','')+"\n")
        for nl in range(0,len(lev2),1):
            perc[nm,nl,:] = 100*var[nm,nl,:]/var[nm,nl,0]
            if nl in [0,1,2]:
                perc2[nm,nl,:] = 100*var[nm,nl,:]/np.sum(var[nm,0:3,:],axis=0)
            elif nl in range(3,7,1):
                perc2[nm,nl,:] = 100*var[nm,nl,:]/np.sum(var[nm,3:7,:],axis=0)
            else:
                perc2[nm,nl,:] = 100*var[nm,nl,:]/np.sum(var[nm,7:10,:],axis=0)

            output = str(var[nm,nl,0])
            for nr in range(1,len(lats),1):
                output = "%s, %d (%.1f%%) (%.1f%%)"%(output,var[nm,nl,nr],perc[nm,nl,nr],perc2[nm,nl,nr])
            f.write("* "+lev2[nl]+", "+output+"\n")
    f.close()

    ds = xr.Dataset(
            {
                "num"  : (["month", "lev", "behv"], var),
                "perc" : (["month", "lev", "behv"], perc),
                "perc2": (["month", "lev", "behv"], perc2),
                },
            coords={
                "month": range(0,len(month),1), 
                "lev"  : (["lev"], lev2),
                "behv" : (["behv"],behv),
                },
            )
    ds.attrs["description"] = "number and percent of cycle, minimum lifetime %d hour"%time
    ds.to_netcdf(fileout,"w")

#===============================================
# draw figure 
#===================================================
if drawannual == 1:
    os.chdir("/home/users/qd201969/uor_track/fig")
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
    fig.savefig(figdir+".png")
    #fig.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.01)

#===============================================
# draw trajectory and filter box 
#===================================================
if drawbox == 1:
    os.chdir("/home/users/qd201969/uor_track/fig")
    if calcbehv != 1:
        fvar = xr.open_dataset(fileout)
        var  = fvar['num'].data
        perc = fvar['perc'].data
        perc2= fvar['perc2'].data

    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    ff = open("%smatch_%s/number"%(path,suffix),'r')
    lines = ff.readlines()
    
    for nr in range(4,len(lats),1):
        if nr == 0:
            suffix2 = ""
        else:
            suffix2 = "_"+behv[nr]
        
        cfp.setvars(file="traj-match2-"+prefix+"_"+suffix+suffix2+".png")
        cfp.gopen(figsize=[20, 20],rows=4,columns=3,wspace=0.1,hspace=0.015,bottom=0.5)
        cfp.mapset(lonmin=0, lonmax=150, latmin=15, latmax=70)
        for nl in range(len(lev2)-1,-1,-1):#,len(f),1):
            np=npp[nl]
            term = lines[nl].strip().split(" ")
            filname = term[0].strip()+suffix2
            print(filname)
            
            cfp.gpos(np)
            cfp.levs(manual=[1500,3000,4500])
            cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2,\
                    title="%s %s %d (%.1f%%) (%.1f%%)"%(lev2[nl],behv[nr],var[0,nl,nr],\
                    perc[0,nl,nr],perc2[0,nl,nr])) #suffix+" "+
            cfp.levs()
            
            for nrr in range(0,len(lats),1):
                cfp.plotvars.mymap.plot([lonl[nrr],lonl[nrr],lonr[nrr],lonr[nrr],lonl[nrr]],
                        [latn[nrr],lats[nrr],lats[nrr],latn[nrr],latn[nrr]], 
                        linewidth=4, color='k', transform=ccrs.PlateCarree()) # filter box

            if var[0,nl,nr] == 0:
                continue
            
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
            
        cfp.cbar(position=[0.2, 0.48, 0.6, 0.01], title='Relative Vorticity (Hz)*1e5')
        cfp.gclose()
        subprocess.run("mogrify -bordercolor white -trim ./traj-match2-"+prefix+"_"+suffix+suffix2+".png",shell=True) 
        #subprocess.run("rm /home/users/qd201969/ERA5-1HR-lev/*.nc",shell=True) 

