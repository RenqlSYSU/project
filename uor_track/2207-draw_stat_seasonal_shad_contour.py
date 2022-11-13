#!/usr/bin/env python
import sys
import subprocess
import xarray as xr
import numpy as np
import gc #garbage collector
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import cmaps
#matplotlib.use('Agg')

title_font=14
label_font=10
plt.rcParams["font.weight"] = "bold"
font = {'family': 'sans-serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        }

lev = [850,500,250]
prefix = "ff"
suffix = 'local'#'remote'# 
#suffix = ''#'_6local'#,'_total']'_6outside'#,
titls= ['DJF','MAM','JJA','SON']
numod= [chr(i) for i in range(97,115)]
figdir = "/home/ys17-23/Extension2/renql/project/uor_track/fig"
path = '/home/ys17-23/Extension2/renql/ERA5-1HR-lev/statistic'
uwndpath = '/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_u_1979-2020.nc'
if suffix in ['_6local','local','remote']:
    lonl=50  #0  #
    lonr=150#360#
    lats=15 #20 #
    latn=55 #90 #
    suf1='local'
if suffix == '_6outside':
    lonl=20  #0  #
    lonr=140#360#
    lats=15 #20 #
    latn=60 #90 #
    suf1='remote'
if suffix == '':
    lonl=0 
    lonr=150
    lats=15 
    latn=70 
    suf1='total'
lat_sp = 20
lon_sp = 30 #60 #

files = '%s/ff_250_1980-2020_stat.nc'%path
f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset(uwndpath)
da = ds['u'].sel(level=200,longitude=ilon,
    latitude=ilat,method="nearest").load()
uwnd = da.groupby(da.time.dt.month).mean('time').data
del ds, da

ds = xr.open_dataset("/home/ys17-23/Extension2/renql/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").data
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

def main_run():
    '''
    # 6local & 6outside
    cnlev1 = np.hstack((np.arange(1,4,1),np.arange(7,20,3)))
    cnlev2 = np.hstack((np.arange(1,8.5,2.5),np.arange(13.5,50,5)))
    draw_shad_cont_seasonal_4x3(suffix,'msp',np.arange(15,85,5),'Speed',
            'gden',cnlev1,'Genesis','lden',cnlev1,'Lysis',False,[4,60])
    draw_shad_cont_seasonal_4x3(suffix,'mstr',np.arange(0,14,1),'Intensity',
            'tden',cnlev2,'Track','lden',cnlev1,'U200',True,[8.5,60])
    '''
    # local match remote or local
    cnlev1 = np.hstack((np.arange(0.5, 2 ,0.5),np.arange( 2 ,20 , 1))) # genesis & lysis
    cnlev2 = np.hstack((np.arange(0.5,3.5, 1 ),np.arange(3.5,50 , 2))) # track
    draw_shad_cont_seasonal_4x3(suffix,'msp',np.arange(15,85,5),'Speed',
            'gden',cnlev1,'Genesis','lden',cnlev1,'Lysis',False,[2,60])
    draw_shad_cont_seasonal_4x3(suffix,'mstr',np.arange(0,14,1),'Intensity',
            'tden',cnlev2,'Track','lden',cnlev1,'U200',True,[3.5,60])

    '''
    # total cyclones
    cnlev1 = np.hstack((np.arange(2,5,1),np.arange(7,20,3)))
    cnlev2 = np.hstack((np.arange(5,20,5),np.arange(26,80,6)))
    draw_shad_cont_seasonal_4x3(suffix,'gden',np.arange(0,7,0.5),'Genesis',
            'tden',cnlev2,'Track','lden',cnlev1,'U200',True,[20,80])
    draw_shad_cont_seasonal_4x3(suffix,'tden',np.arange(0,42,3),'Track',
            'gden',cnlev1,'Gensis','lden',cnlev1,'U200',True,[5,80])
    '''

def draw_shad_cont_seasonal_4x3(suffix,varname,cnlev,label,
    varname1,cnlev1,label1,varname2,cnlev2,label2,jetoption,dash):
    
    nrow = 4 #6 #
    ncol = 3 #2 #
    bmlo = 0.4 #0.25 #

    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
    
    colr = cmaps.topo_15lev(range(0,16,1))[::-1]
    print(colr)
    colr[8,:] = colors.to_rgba('y')
    #colr[8:13,:] = colr[9:14,:]
    colr[7,:] = colors.to_rgba('c')
    ncmap = colors.ListedColormap(colr)
    norm = colors.BoundaryNorm(boundaries=cnlev,
        ncolors=ncmap.N,extend='both')
    
    for nl in range(0,len(lev),1):
        files = '/home/ys17-23/Extension2/renql/project/uor_track/mdata/statistic/ff_match_%dlocal%s_6dist.nc'%(lev[nl],suffix)
        #files = '%s/%s_%d_1980-2020%s_stat.nc'%(path,prefix,lev[nl],suffix)
        f = xr.open_dataset(files)
        print("")
        print(files)
        var = read_stat(f,varname,lev[nl])
        var1 = read_stat(f,varname1,lev[nl])
        var2 = read_stat(f,varname2,lev[nl])
    
        for nm in range(0,nrow,1):
            if nm == 0:
                shad = (var[0,:,:]+var[1,:,:]+var[11,:,:])/3.0
                cont1 = (var1[0,:,:]+var1[1,:,:]+var1[11,:,:])/3.0
                cont2 = (var2[0,:,:]+var2[1,:,:]+var2[11,:,:])/3.0
                uwnd1 = (uwnd[0,:,:]+uwnd[1,:,:]+uwnd[11,:,:])/3.0
            else:
                shad = np.mean(var[(3*nm-1):(3*nm+2),:,:],axis=0)
                cont1 = np.mean(var1[(3*nm-1):(3*nm+2),:,:],axis=0)
                cont2 = np.mean(var2[(3*nm-1):(3*nm+2),:,:],axis=0)
                uwnd1 = np.mean(uwnd[(3*nm-1):(3*nm+2),:,:],axis=0)
            axe = ax[nm][nl]
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
                edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title("(%s) %dhPa %s %s"%(numod[3*nm+nl],lev[nl],titls[nm],
                suffix.strip('_')),fontsize=title_font,fontdict=font)

            cont = axe.contourf(ilon, ilat, shad, cnlev, 
                 transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000,4500], linestyles='solid', 
                 transform=ccrs.PlateCarree(),colors='black',linewidths=2.0)

            line1 = axe.contour(ilon, ilat, cont1, cnlev1, linestyles='solid',
                 transform=ccrs.PlateCarree(),colors='b',linewidths=1.5)
            if not(varname1=='gden' and nl==2 and nm==3):
                line3 = axe.contour(ilon, ilat, cont1, dash, linestyles='dashed', 
                     transform=ccrs.PlateCarree(),colors='b',linewidths=1.5)
            
            if jetoption :
                jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], linestyles='solid', 
                     transform=ccrs.PlateCarree(),colors='r',linewidths=2.0)
            else :
                line2 = axe.contour(ilon, ilat, cont2, cnlev2, linestyles='solid', 
                     transform=ccrs.PlateCarree(),colors='r',linewidths=1.5)
                line4 = axe.contour(ilon, ilat, cont2, dash, linestyles='dashed', 
                     transform=ccrs.PlateCarree(),colors='r',linewidths=1.5)
            
            if nl == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nm == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

#    legend_elements = [Line2D([0], [0], color='r', lw=4, label=label1),
#         Line2D([0], [0], color='m', lw=4, label=label2)]
#    ax[3][0].legend(handles=legend_elements, bbox_to_anchor=(0.1, bmlo+0.005, 0.15, 0.01), 
#        bbox_transform=fig.transFigure, ncol=2)

    position = fig.add_axes([0.45, bmlo+0.005, 0.45, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.figtext(0.35,bmlo-0.005, label,fontsize=title_font,
        horizontalalignment='left',verticalalignment='bottom')

    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig('%s/stat_seasonal%s_%s.png'%(figdir,suffix,varname), 
        bbox_inches='tight',pad_inches=0.01)

def draw_shad_cont_3x1(suffix,varname,cnlev,label,
    varname1,cnlev1,label1,varname2,cnlev2,label2,jetoption,dash):
    
    nrow = 3 #6 #
    ncol = 1 #2 #
    bmlo = 0.5

    fig = plt.figure(figsize=(12,12),dpi=150)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(
        projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True

    ncmap = colors.ListedColormap(cmaps.topo_15lev(range(0,16,1))[::-1])
    norm = colors.BoundaryNorm(boundaries=cnlev,
        ncolors=ncmap.N,extend='both')

    uwnd1 = uwnd.mean(axis=0)
    for nr in range(0,len(lev),1):
        files = '%s/%s_%d_1980-2020%s_stat.nc'%(path,prefix,lev[nr],suffix)
        f = xr.open_dataset(files)
        print("")
        print(files)
        shad = np.nanmean(read_stat(f,varname,lev[nr]),axis=0)
        cont1 = np.nanmean(read_stat(f,varname1,lev[nr]),axis=0)
        cont2 = np.nanmean(read_stat(f,varname2,lev[nr]),axis=0)
    
        axe = ax[nr]
        axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],
            edgecolor='k'), linewidth=0.8, zorder=1)
        axe.set_title("(%s) %dhPa %s"%(numod[nr],lev[nr],label
            ),fontsize=title_font,fontdict=font)

        cont = axe.contourf(ilon, ilat, shad, cnlev, 
             transform=ccrs.PlateCarree(),cmap=ncmap,extend='both',norm=norm)
        topo = axe.contour(ilon, ilat, phis, [1500,3000], linestyles='solid', 
             transform=ccrs.PlateCarree(),colors='black',linewidths=0.8)

        line1 = axe.contour(ilon, ilat, cont1, cnlev1, linestyles='solid',
             transform=ccrs.PlateCarree(),colors='b',linewidths=1.2)
        line3 = axe.contour(ilon, ilat, cont1, dash, linestyles='dashed', 
             transform=ccrs.PlateCarree(),colors='b',linewidths=1.5)
    
        if jetoption :
            jets = axe.contour(ilon, ilat, uwnd1, [30,40,50], linestyles='solid', 
                 transform=ccrs.PlateCarree(),colors='r',linewidths=2.0)
        else :
            line2 = axe.contour(ilon, ilat, cont2, cnlev2, linestyles='solid', 
                 transform=ccrs.PlateCarree(),colors='r',linewidths=1.2)
            line4 = axe.contour(ilon, ilat, cont2, dash, linestyles='dashed', 
                 transform=ccrs.PlateCarree(),colors='r',linewidths=1.5)
        
        axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
        axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
        if nr == (nrow-1):
            axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
            axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))
    
    position = fig.add_axes([0.35, bmlo+0.001, 0.3, 0.01]) #left, bottom, width, height
    cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig('%s/stat_seasonal%s_%s.png'%(figdir,suffix,varname), 
        bbox_inches='tight',pad_inches=0.01)

def read_stat(f,varname,nl):
    var = f[varname].sel(long=ilon,lat=ilat).load()
    if varname == 'mten':
        var.data = var.data*24
    if varname[-3:] != 'den':
        tden = f['tden'].sel(long=ilon,lat=ilat).data
        mask = tden < 1.0
        var.data = np.ma.array(var.data,mask=mask)
        del mask, tden
    
    if nl==850:
        var.data = np.ma.array(var.data,mask=(
            np.broadcast_to(phis, var.shape)>1500))
    
    print('%s: min:%f ; max:%f'%(var.long_name,
        np.nanmin(var.data),np.nanmax(var.data)))
    return var.data
    
if __name__=='__main__':
    main_run()

