#!/usr/bin/env python
'''
read uwnd, vwnd, z to draw monthly wind (vector) and geopotential (shaded)
Loop through the height (850, 500, 250)
maybe later the shaded variable can be changed for t, PV, dtdy

20211007
'''
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
from scipy import stats

def draw_vect_shad(month,var,uwnd,vwnd,pvar,pu,pv):
    siglvl = 0.1
    usig = uwnd
    vsig = vwnd
    mask1 = np.array([pu>siglvl,pv>siglvl]).all(axis=0)
    usig.values = np.ma.array(usig.values,mask=mask1)
    vsig.values = np.ma.array(vsig.values,mask=mask1)
    mask2 = np.array([pu<=siglvl,pv<=siglvl]).any(axis=0)
    uwnd.values = np.ma.array(uwnd.values,mask=mask2)
    vwnd.values = np.ma.array(vwnd.values,mask=mask2)
    print(usig)
    print(uwnd)

    xbar=[0.04,0.36,0.68]
    for nl1 in range(0,3,1):
        fig = plt.figure(figsize=(12,12),dpi=300)
        ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
        for nr in range(0,nrow,1):
            nb = nr+1
            for nc in range(0,ncol,1):
                nl = nc
                if cnlvl[nl][0] < 0 :
                    fcolors = cmaps.BlueDarkRed18
                else:
                    fcolors = cmaps.precip2_17lev
                cnlevels = np.arange(cnlvl[nl][0], cnlvl[nl][0]+cnlvl[nl][1]*(fcolors.N-1), cnlvl[nl][1])
                norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')

                axe = ax[nr][nc]
                axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
                axe.set_title("%s %d TE %s u_v_z(%dhPa)"%(month,lev[nl1],behv[nb],levc[nl]),fontsize=SMFONT)

                shad = axe.contourf(ilon, ilat, var[nl1,nb,nl,:,:], cnlevels,
                        transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
                sigf = axe.contourf(ilon, ilat, pvar[nl1,nb,nl,:,:],[siglvl,1], 
                        hatches=['xxx',None,None],colors="none",extend='both',transform=ccrs.PlateCarree())
                
                wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd[nl1,nb,nl,::q_mis,::q_mis],vwnd[nl1,nb,nl,::q_mis,::q_mis],
                        pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="gray",zorder=4,
                        width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())
                #wsig = axe.quiver(ilon[::q_mis], ilat[::q_mis], usig[nl1,nb,nl,::q_mis,::q_mis],vsig[nl1,nb,nl,::q_mis,::q_mis],
                #        pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="black",zorder=4,
                #        width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

                topo = axe.contour(ilon, ilat, phis, [1500,3000],
                             transform=ccrs.PlateCarree(),colors='black',linewidths=1.2)

                if dbox >= 1 :
                    axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                             linewidth=2, color='black', transform=ccrs.PlateCarree()) # filter box

                if nc == 0:
                    axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                    axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
                if nr == (nrow-1):
                    axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                    axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
                    axe.quiverkey(wind, 0.95, -0.3, vcref[nl], r'$%d m/s$'%vcref[nl], labelpos='N',coordinates='axes')
                    position = fig.add_axes([xbar[nc], bmlo, 0.28, 0.008]) #left, bottom, width, height
                    cb = plt.colorbar(shad, cax=position ,orientation='horizontal')#, shrink=.9)

        #plt.figtext(0.02,bmlo-0.005, "%dhPa TE"%(lev[nl]), fontsize=MIDFONT,
        #        horizontalalignment='left',verticalalignment='bottom')
        plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
        plt.savefig(figdir+"regr_interannual_%s_%d_%s.png"%(drawvar[0],lev[nl1],month), bbox_inches='tight',pad_inches=0.01)

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30

nrow = 4
ncol = 3
bmlo = 0.4
BIGFONT=22
MIDFONT=14
SMFONT=10

lev = [850,500,250]
levc = [850,500,250]
behv = ["ALL" ,"NTN" ,"STN" ,"PAS" ,"LYS" ]#,"DIF"]
filname = '/home/users/qd201969/uor_track/mdata/regr_interannual_' # t.nc'
varname = ['regr','prob']
drawvar = ['z']
unit    = ['m']
vcref =[2,4,5] # different levels 
cnlvl =[[-40 ,5  ],
        [-120,15 ],
        [-200,25]]
q_mis=5
dbox = 1
flats = 27 #int(sys.argv[2])
flatn = 45 #int(sys.argv[3])
flonl = 60 #int(sys.argv[4])
flonr = 90 #int(sys.argv[5])
figdir = "/home/users/qd201969/uor_track/fig/"
months = ["DJF","MAM","JJA","SON"]

f = xr.open_dataset(filname+'z.nc')
lat = f.lat
lon = f.lon
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
var  = f[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
pvar = f[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()
var.data = var.data/9.8

ds = xr.open_dataset(filname+'u.nc')
uwnd = ds[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
pu   = ds[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()

ds = xr.open_dataset(filname+'v.nc')
vwnd = ds[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
pv   = ds[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

for nm in range(0,len(months),1):
    draw_vect_shad(months[nm],var[:,:,nm,:,:,:],uwnd[:,:,nm,:,:,:],vwnd[:,:,nm,:,:,:],\
            pvar[:,:,nm,:,:,:],pu[:,:,nm,:,:,:],pv[:,:,nm,:,:,:])

