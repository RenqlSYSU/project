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

def draw_vect_shad(month,var,uwnd,vwnd,numb,vari,uvari,vvari):
    xbar=[0.04,0.36,0.68]
    siglvl = 0.1
    for nl1 in range(0,3,1):
        fig = plt.figure(figsize=(12,12),dpi=300)
        ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree())) #sharex=True, sharey=True
        for nr in range(0,nrow,1):
            nb = nr+1
            if diff == 1:
                #t,pvar = stats.ttest_ind_from_stats(var[nl1,nb,:,:,:],vari[nl1,nb,:,:,:],numb[nl1,nb],
                #        var[nl1,0,:,:,:],vari[nl1,0,:,:,:],numb[nl1,0],equal_var=True)
                #t,pu = stats.ttest_ind_from_stats(uwnd[nl1,nb,:,:,:],uvari[nl1,nb,:,:,:],numb[nl1,nb],
                #        uwnd[nl1,0,:,:,:],uvari[nl1,0,:,:,:],numb[nl1,0],equal_var=True)
                #t,pv = stats.ttest_ind_from_stats(vwnd[nl1,nb,:,:,:],vvari[nl1,nb,:,:,:],numb[nl1,nb],
                #        vwnd[nl1,0,:,:,:],vvari[nl1,0,:,:,:],numb[nl1,0],equal_var=True)
                #print(pvar)
                
                var[nl1,nb,:,:,:] = (var[nl1,nb,:,:,:] - var[nl1,0,:,:,:])/9.8
                uwnd[nl1,nb,:,:,:]=uwnd[nl1,nb,:,:,:] - uwnd[nl1,0,:,:,:]
                vwnd[nl1,nb,:,:,:]=vwnd[nl1,nb,:,:,:] - vwnd[nl1,0,:,:,:]
                
                #var[nl1,nb,:,:,:].values=np.ma.array(var[nl1,nb,:,:,:].values,mask=(pvar>siglvl))
                
                #mask = np.array([pu>siglvl,pv>siglvl]).all(axis=0)
                #uwnd[nl1,nb,:,:,:].values=np.ma.array(uwnd[nl1,nb,:,:,:].values,mask=mask)
                #vwnd[nl1,nb,:,:,:].values=np.ma.array(vwnd[nl1,nb,:,:,:].values,mask=mask)
            
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
                axe.set_title("%s %d TE %s (%d) u_v_z(%dhPa)"%(month,lev[nl1],behv[nb],numb[nl1,nb],levc[nl]),fontsize=SMFONT)

                shad = axe.contourf(ilon, ilat, var[nl1,nb,nl,:,:], cnlevels,
                             transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
                
                wind = axe.quiver(ilon[::q_mis], ilat[::q_mis], uwnd[nl1,nb,nl,::q_mis,::q_mis],vwnd[nl1,nb,nl,::q_mis,::q_mis],
                        pivot='mid',units='inches',scale=vcref[nl]*3,scale_units='inches',color="dimgray",
                        width=0.02,headwidth=3,headlength=4.5,transform=ccrs.PlateCarree())

                #cont = axe.contour(ilon, ilat, var[nm,:,:], np.arange(1000,15000,cnlvl[nl][1]), 
                #             transform=ccrs.PlateCarree(), colors='darkviolet', linewidths=2)

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
        plt.savefig(figdir+"composite_%s_%d_%s_diff%d.png"%(drawvar[0],lev[nl1],month,diff), bbox_inches='tight',pad_inches=0.01)

lonl=0  #0  #
lonr=150#360#
lats=15 #0  #
latn=70 #90 #
lat_sp = 20
lon_sp = 30

diff=1 # if diff=1, then draw diff
nrow = 6
ncol = 3
bmlo = 0.25
BIGFONT=22
MIDFONT=14
SMFONT=10

lev = [850,500,250]
levc = [850,500,250]
behv = ["PAS" ,"NTP" ,"STP" ,"NTL" ,"STL" ,"LYS" ]#,"DIF"]
#filname = '/home/users/qd201969/uor_track/mdata/regr_interannual_' # t.nc'
#varname = ['regr','prob']
filname = '/home/users/qd201969/uor_track/mdata/comp_6h_season_daily00_' # t.nc'
varname = ['var','vari']
drawvar = ['z']
unit    = ['m']
if diff == 1:
    vcref =[5,8,10] # different levels 
    cnlvl =[[-40 ,5  ],
            [-120,15 ],
            [-200,25]]
else:
    vcref =[10,20,30] # different levels 
    cnlvl =[[1300 ,20 ],
            [5050 ,60 ],
            [9300 ,100]]
q_mis=5
dbox = 1
flats = 25 #int(sys.argv[2])
flatn = 45 #int(sys.argv[3])
flonl = 60 #int(sys.argv[4])
flonr = 105 #int(sys.argv[5])
figdir = "/home/users/qd201969/uor_track/fig/"
months = ["DJF","MAM","JJA","SON"]

f = xr.open_dataset(filname+'z.nc')
lat = f.lat
lon = f.lon
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]
numb = f['numb']
var = f[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
vari = f[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()
var.data = var.data
print("var[0,0,1,2,:,:]")
print(var[0,0,1,2,:,:])
print("var[0,1,1,2,:,:]")
print(var[0,1,1,2,:,:])

ds = xr.open_dataset(filname+'u.nc')
uwnd = ds[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
uvari = ds[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()

ds = xr.open_dataset(filname+'v.nc')
vwnd = ds[varname[0]].sel(level=levc,lon=ilon,lat=ilat).load()
vvari = ds[varname[1]].sel(level=levc,lon=ilon,lat=ilat).load()
print("uwnd")
print(uwnd[0,0,1,2,:,:])
print(uwnd[0,1,1,2,:,:])
print("vwnd")
print(vwnd[0,0,1,2,:,:])
print(vwnd[0,1,1,2,:,:])

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

for nm in range(0,len(months),1):
    draw_vect_shad(months[nm],var[:,:,nm,:,:,:],uwnd[:,:,nm,:,:,:],vwnd[:,:,nm,:,:,:],numb[:,:,nm],\
            vari[:,:,nm,:,:,:],uvari[:,:,nm,:,:,:],vvari[:,:,nm,:,:,:])

