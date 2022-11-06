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
import matplotlib.ticker as mtick
import cartopy.crs as ccrs
from renql import cyc_filter, monthly_calc, life_intensity

title_font=18
label_font=14
plt.rcParams["font.weight"] = "bold"
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black',
        }

if len(sys.argv) < 2 :
    option=2 #int(sys.argv[1]) #Genesis (0)/Lysis (1)/Passing(2)/Passing Time(3)/All Times(4)
    flats = 45  #int(sys.argv[2])
    flatn = 45  #int(sys.argv[3])
    flonl = 60  #int(sys.argv[4])
    flonr = 110 #int(sys.argv[5])
    time = 24 # threshold, hour
    prefix = "fftadd"
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

suffix=str(option)+"_"+str(flats)+str(flatn)+"-"+str(flonl)+str(flonr)
figdir = "/home/users/qd201969/uor_track/fig/"
fileout="/home/users/qd201969/uor_track/mdata/behv4_month_%dh_%s.nc"%(time,suffix)

flonr2 = 105
flatn2 = 45
flats2 = 25
behv = ["ALL" ,"EPAS","ELYS","WPAS","WNTP","WNTL","WLYS" ]#,"DIF"]
lats = [flats ,flats ,flats ,flats2,flatn2,flatn ,flats2]
latn = [flatn ,flatn ,flatn ,flatn2,flatn2,flatn ,flatn2]
lonl = [flonl ,85    ,85    ,flonr2,flonr2,flonl ,flonl ]
lonr = [flonr ,flonr ,flonr ,flonr2,flonr2,flonr2,flonr2]
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

def main_run():
    #calcbehv()
    #drawannual()
    #draw_table()
    draw_3var_distri()
    #drawbox()
    #draw_hist()

def calcbehv():
    var  = np.empty( [len(month),len(lev),len(behv)],dtype=int )  
    perc = np.empty( [len(month),len(lev),len(behv)],dtype=float )  
    for nl in range(0,len(lev),1):
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0+"_"+suffix#+"_"+str(time)

        if not os.path.isfile(filname) :
            filname0 = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0
            var[0,nl,0] = cyc_filter.line_filt(filname0,flats,flatn,flonl,flonr,time,1,"north")

        var[0,nl,:] = cyc_filter.north_behavior(filname,
                lats[1:len(behv)],latn[1:len(behv)],lonl[1:len(behv)],lonr[1:len(behv)],lonx=105)

        for nr in range(0,len(behv),1):
            if nr == 0:
                filname1 = filname
            else:
                filname1 = filname+"_"+behv[nr]
            #var[1:len(month),nl,nr] = monthly_calc.calc_month(filname1,time)
            #var[1:len(month),nl,nr] = monthly_calc.draw_month_traj(
            #   frae, month[1:len(month)],nday[1:len(month)],filname1,lats,latn,lonl,lonr)
    var[0,:,:] = np.sum(var[1:len(month),:,:],axis=0)
    '''
    f = open("%sbehv3_month_%dh_%s"%(figdir,time,suffix),'w')
    f.write(path+prefix+"_1980-2020_"+suffix+"\n")
    f.write("*"*50+"\n")
    f.write("lats"+str(lats)+"\n")
    f.write("latn"+str(latn)+"\n")
    f.write("lonl"+str(lonl)+"\n")
    f.write("lonr"+str(lonr)+"\n")
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
    '''
def draw_hist():
    nrow = len(lev)
    ncol = len(behv)-1
    bmlo = 0.1
    
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
                ax[nl,nr-1].set_xlabel("lifetime (days)",fontsize=label_font,fontdict=font)

            if nr == 1:
                ax[nl,nr-1].set_ylabel("Mean intensity ($10^{-5} s^{-1}$)",fontsize=label_font,fontdict=font)

    position = fig.add_axes([0.1, bmlo+0.05, 0.8, 0.02]) #left, bottom, width, height
    cb = plt.colorbar(patches, cax=position ,orientation='horizontal')#, shrink=.9)
    fig.tight_layout(rect=(0,bmlo,0.94,1)) #w_pad=0.5,h_pad=0.001) #,
    fig.savefig("%slife_inte_%s.png"%(figdir,suffix), bbox_inches='tight',pad_inches=0.08)

def draw_3var_distri():
    #varname = ["Lifetime(day)","Intensity ($10^{-5} s^{-1}$)","Distance(km)"]
    #varname = ["close1 (%)","close2 (%)","GeopoHeight (m)"]
    varname = ["10mWind (m/s)","MaxRain (mm/d)","GeopoHeight (m)"]
    nrow = len(lev)
    ncol = 3
    bmlo = 0.4
    #xbin = [np.arange(0,101,5), # num=(end-start)/inv
    #        np.arange(0,101,5), # close1,close2,geopotion  
    #        np.arange(0,8001,400)]
    xbin = [np.arange(2,22.1,1), # num=(end-start)/inv
            np.arange(0,50.1,2.5), # 10mwind, maxrain, raintime 
            np.arange(0,101,5)]
    #xbin = [np.arange(1,21, 1 ),
    #        np.arange(1,21, 1 ),
    #        np.arange(0,10000,500)]
    var = np.empty( [len(varname),len(behv),len(xbin[0])-1],dtype=float )   
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for nl in range(0,len(lev),1):
        if nl==0 : # for geopotential
            xbin[2]=np.arange(1200,1701,25)
        elif nl==1 :
            xbin[2]=np.arange(5000,6001,50)
        elif nl==2 :
            xbin[2]=np.arange(9500,11001,75)
        
        numb = []
        filname  = path+prefix+"_"+str(lev[nl])+"_1980-2020"+suffix0+"_"+suffix#+"_"+str(time)
        for nb in range(1,len(behv),1):
            if nb == 0:
                filname1 = filname
            else:
                filname1 = filname+"_"+behv[nb]
            
            #life1, inte1, dist1, numb1 = life_intensity.calc_life_intensity(
            #        filname1,flats=25,flatn=45,flonl=57,flonr=110)
            life1, inte1, dist1, numb1 = life_intensity.calc_close_cyclone(
                    filname1,flats=20,flatn=90,flonl=50,flonr=130)
            print("%s :%d"%(filname1,numb1))
            var[0,nb,:],term = np.histogram(life1,xbin[0])
            var[1,nb,:],term = np.histogram(inte1,xbin[1])
            var[2,nb,:],term = np.histogram(dist1,xbin[2])
            var[:,nb,:] = var[:,nb,:]*100/numb1
            numb.append(numb1)
                
        for nv in range(0,3,1):
            ax[nl][nv].grid(True, which="both", color='grey', linestyle='--', linewidth=1)
            ax[nl][nv].yaxis.set_major_formatter(mtick.FormatStrFormatter('%i'))
            for nb in range(1,len(behv),1):
                ax[nl][nv].plot(xbin[nv][1:],var[nv,nb,:],linewidth=2)
            
            if nl==2 :
                ax[nl][nv].set_xlabel(varname[nv],fontsize=label_font-4,fontdict=font)
            if nv==0 :
                ax[nl][nv].set_ylabel("%dhPa percent"%lev[nl],fontsize=label_font-4,fontdict=font)

    ax[0][2].legend(behv[1:len(behv)], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%slife_inte_dist_%dh_%s"%(figdir,time,suffix))

def drawannual():
    nrow = 3
    ncol = 2
    bmlo = 0.4
    fvar = xr.open_dataset(fileout)
    var = fvar['num'].data
    perc= fvar['perc'].data
    var = np.rint(var/41)

    fig = plt.figure(figsize=(9,9),dpi=200)
    ax = fig.subplots(nrow, ncol) #sharex=True, sharey=True
    for nl in range(0,len(lev),1):
        ax[nl][0].set_title(str(lev[nl])+" number "+suffix,fontsize=title_font,fontdict=font)
        ax[nl][1].set_title(str(lev[nl])+" percent "+suffix,fontsize=title_font,fontdict=font)
        ax[nl][0].set_xlim(1, 12)
        ax[nl][1].set_xlim(1, 12)
        for nr in range(1,len(lats),1):
            ax[nl][0].plot(imonth,var[1:len(month),nl,nr],linewidth=2)
            ax[nl][1].plot(imonth,perc[1:len(month),nl,nr],linewidth=2)
        ax[nl][0].grid(True, which="both", color='grey', linestyle='--', linewidth=1)
        ax[nl][1].grid(True, which="both", color='grey', linestyle='--', linewidth=1)
    
    ax[2][1].legend(behv[1:len(behv)], loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sbehv3_month_%dh_%s"%(figdir,time,suffix))
    #fig.savefig("%sbehv3_month_%dh_%s"%(figdir,time,suffix), bbox_inches='tight',pad_inches=0.01)
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    ax  = fig.subplots(1, 1) #sharex=True, sharey=True
    ax.set_title("total number "+suffix,fontsize=title_font,fontdict=font)
    for nl in range(0,len(lev),1):
        ax.plot(imonth,var[1:len(month),nl,0],linewidth=2,label=str(lev[nl]))
    
    ax.legend(loc='upper right')
    fig.tight_layout(w_pad=0.5,h_pad=1) #,rect=(0,bmlo,1,1)
    fig.savefig("%sbehv3_month_%dh_%s_total.png"%(figdir,time,suffix))

def drawbox():
    bmlo = 0.25
    fvar = xr.open_dataset(fileout)
    var = fvar['num'].data
    perc= fvar['perc'].data

    f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
    phis=f0[2]
    print(repr(phis))
    phis=phis/9.8 # transfer from m2/s2 to m
    
    cfp.setvars(file="traj-"+prefix+"_"+suffix+".png")
    cfp.gopen(figsize=[20, 20],rows=6,columns=3,wspace=0.1,hspace=0.02,bottom=bmlo) #0.55
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

            if var[0,nl,nr] in [0,1]:
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
            subprocess.run("rm %s.nc"%filname,shell=True) 

            cfp.cscale('precip2_17lev')
            cfp.levs(min=0.0, max=8.0, step=0.5)
            cfp.traj(g, zorder=0, legend_lines=True, colorbar=False, linewidth=1.5)
            
    cfp.cbar(position=[0.2, bmlo-0.02, 0.6, 0.01], title='Relative Vorticity (Hz)*1e5') #0.53
    cfp.gclose()
    subprocess.run("mogrify -bordercolor white -trim ./traj-"+prefix+"_"+suffix+".png",shell=True) 

def draw_table():
    fvar = xr.open_dataset(fileout)
    var = fvar['num'][0,:,:].data
    perc= fvar['perc'][0,:,:].data
    
    col = behv
    row = ['%d hPa' %x for x in lev]
    vals = []
    for nl in range(len(lev)):
        vals.append(['%d (%.2f%%)'%(x,y) for x,y in zip(var[nl,:],perc[nl,:])])
        vals[nl][0] = str(var[nl,0])

    fig = plt.figure(figsize=(9,9),dpi=200)
    tab = plt.table(cellText=vals, colLabels=col, rowLabels=row,
        loc='center', cellLoc='center',rowLoc='center',colLoc='center')
        #rowColours='grey',colColours='grey')
    tab.scale(1,2) 
    plt.axis('off')
    fig.savefig("%sbehv4_table_%s.png"%(figdir,suffix), bbox_inches='tight',pad_inches=0.08)
    subprocess.run("mogrify -bordercolor white -trim %sbehv4_table_%s.png"%(figdir,suffix),shell=True) 

if __name__=='__main__':
    main_run()
