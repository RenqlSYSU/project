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
#matplotlib.use('Agg')

lonl=0  #0  #
lonr=150#360#
lats=15 #20 #
latn=70 #90 #
lat_sp = 20
lon_sp = 30 #60 #
nrow = 4 #6 #
ncol = 3 #2 #
bmlo = 0.4 #0.25 #
title_font=14
label_font=10

# level == 0 #small range
lev=[[0    ,320  ,20  ], # 0Feature Density
     [0    ,3.2  ,0.2 ], # 1Genesis Density
     [0    ,3.2  ,0.2 ], # 2Lysis Density
     [-1   ,1    ,1   ], # 3Mean Area
     [-0.8 ,0.8  ,0.1 ], # 4Mean Growth/Decay Rate
     [-1   ,1    ,1   ], # 5Mean Anisotropy
     [0    ,8    ,0.5 ], # 6Mean Lifetime
     [0    ,80   ,5   ], # 7Mean Speed
     [0    ,16   ,1   ], # 8Mean Intensity
     [-1.6 ,1.6  ,0.2 ], # 9Mean Tendency
     [-1   ,1    ,1   ], # 10Spare1
     [-1   ,1    ,1   ], # 11Spare2
     [0    ,80   ,5   ], # 12Std of Speed
     [0    ,3.2  ,0.2 ], # 13Std of Intensity
     [0    ,16   ,1   ], # 14Track Density
     [-1   ,1    ,1   ], # 15X-component of Mean Orientation Vector
     [-40  ,40   ,5   ], # 16X-component of Mean Velocity
     [-1   ,1    ,1   ], # 17Y-component of Mean Orientation Vector
     [-40  ,40   ,5  ]] # 18Y-component of Mean Velocity

draw_var = ["fden","gden","lden","marea","mgdr","",
            "mlif","msp" ,"mstr","mten" ,""    ,"",
            ""    ,""    ,"tden",""    ] # 7 variables
#draw=[8,9,6]
draw=[1,2,14]
#draw=[14]
#draw=[1,2,14,8,9,6]

if len(sys.argv) < 2 :
    filename = 'ff_250_1980-2020' #_2_2545-6080
    files = '/home/users/qd201969/ERA5-1HR-lev/statistic/'+filename+'_stat.nc'
    level = 2
    figtitle = '250'#_2_2545-6080
    dbox = 0 
else:
    filename = sys.argv[1]  #'ff_250_500_no'
    files = sys.argv[2] #'/home/users/qd201969/ERA5-1HR-lev/match'+filt+'/statistic/'+filename+'_stat_'
    level = int(sys.argv[3])
    figtitle = sys.argv[4]
    dbox = int(sys.argv[5])

if dbox >= 1 :
    flats = int(sys.argv[6])
    flatn = int(sys.argv[7])
    flonl = int(sys.argv[8])
    flonr = int(sys.argv[9])

if level == 1: #middle range 
    lev[ 1][:]=[0    ,4.8  ,0.3 ]
    lev[ 2][:]=[0    ,4.8  ,0.3 ]
    lev[14][:]=[0    ,24   ,1.5 ]
if level == 2: # large range,use total level bar 
    lev[ 0][:]=[0    ,1600 ,100 ]
    lev[ 1][:]=[0    ,8    ,0.5 ]
    lev[ 2][:]=[0    ,8    ,0.5 ]
    lev[14][:]=[0    ,48   ,3   ]
titls=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

f = xr.open_dataset(files)
lat = f.lat
lon = f.long
ilon = lon[(lon>=lonl) & (lon<=lonr)]
ilat = lat[(lat>=lats) & (lat<=latn)]

ds = xr.open_dataset('/home/users/qd201969/data/ERA5_mon_u_1979-2020.nc')
da = ds['u'].sel(level=200,longitude=ilon,latitude=ilat,method="nearest").load()
# increased performance by loading data into memory first, e.g., with load()
uwnd = da.groupby(da.time.dt.month).mean('time')
del ds, da

ds = xr.open_dataset("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis = ds['PHIS'].sel(lon=ilon,lat=ilat,method="nearest").load()
phis = phis/9.8 # transfer from m2/s2 to m
del ds
gc.collect()

figdir = "/home/users/qd201969/uor_track/fig/month_"+filename
for nv in range(0,len(draw),1):#,len(f),1):
    var = f[draw_var[draw[nv]]].sel(long=ilon,lat=ilat).load()
    print(var)
    if draw[nv] == 9:
        var=var*24
    
    if draw[nv] > 2 and draw[nv] != 14:
        tden = f['tden'].sel(long=ilon,lat=ilat).load()
        mask = tden < 1.0
        var.values=np.ma.array(var.values,mask=mask)
    
    cnlevels = np.arange(lev[draw[nv]][0], lev[draw[nv]][1]+lev[draw[nv]][2], lev[draw[nv]][2])
    if lev[draw[nv]][0] < 0 :
        fcolors = cmaps.BlueDarkRed18
    else:
        fcolors = cmaps.precip2_17lev
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(12,12),dpi=300)
    ax = fig.subplots(nrow, ncol, subplot_kw=dict(projection=ccrs.PlateCarree(central_longitude=180.0))) #sharex=True, sharey=True
    nm = -1
    #for nm in range(0,len(titls),1):
    for nr in range(0,nrow,1):
        for nc in range(0,ncol,1):
            nm = nm+1 
            axe = ax[nr][nc]
            #axe = plt.subplot(4,3,nm+1,projection=ccrs.PlateCarree())    #创建子图
            #axe.add_feature(cfeat.COASTLINE.with_scale('110m'),edgecolor='black', linewidth=0.8, zorder=1) 
            axe.add_feature(cfeat.GSHHSFeature(levels=[1,2],edgecolor='k'), linewidth=0.8, zorder=1)
            axe.set_title(figtitle+" "+titls[nm],fontsize=title_font)

            cont = axe.contourf(ilon, ilat, var[nm,:,:], cnlevels, 
                         transform=ccrs.PlateCarree(),cmap=fcolors,extend='both',norm=norm)
            topo = axe.contour(ilon, ilat, phis, [1500,3000], 
                         transform=ccrs.PlateCarree(),colors='black',linewidths=1.5)
            if dbox >= 1 :
                axe.plot([flonl,flonl,flonr,flonr,flonl],[flatn,flats,flats,flatn,flatn], 
                         linewidth=2.5, color='black', transform=ccrs.PlateCarree()) # filter box
            jets = axe.contour(ilon, ilat, uwnd[nm,:,:], [30,40,50], 
                         transform=ccrs.PlateCarree(),colors='darkviolet',linewidths=2)
            if nc == 0:
                axe.set_yticks(np.arange(lats,latn,lat_sp), crs=ccrs.PlateCarree())
                axe.yaxis.set_major_formatter(LatitudeFormatter(degree_symbol=''))
            if nr == (nrow-1):
                axe.set_xticks(np.arange(lonl,lonr,lon_sp), crs=ccrs.PlateCarree())
                axe.xaxis.set_major_formatter(LongitudeFormatter(degree_symbol=''))

    #fig.subplots_adjust(left=0.2,bottom=bmlo) # wspace control horizontal space
    if (nrow/ncol) >= 2.0: 
        position = fig.add_axes([0.99, bmlo+0.05, 0.01, 0.6]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='vertical')#, shrink=.9)
        cb.set_label(label=var.long_name, size=title_font) #, weight='bold'
    else:
        position = fig.add_axes([0.2, bmlo+0.005, 0.7, 0.01]) #left, bottom, width, height
        cb = plt.colorbar(cont, cax=position ,orientation='horizontal')#, shrink=.9)
        plt.figtext(0.02,bmlo-0.005, var.long_name,fontsize=title_font,
                horizontalalignment='left',verticalalignment='bottom')
    #axe.text(0.5,1.2, filename+' '+var.long_name,fontsize=title_font,
    #        horizontalalignment='center',verticalalignment='bottom',transform=ax[0][1].transAxes)
    #plt.suptitle(filename+' '+var.long_name,x=0.5,y=0.9,fontsize=title_font/2)
    plt.tight_layout(w_pad=0.5,rect=(0,bmlo,1,1))
    plt.savefig(figdir+var.long_name.replace(" ","")+".png", bbox_inches='tight',pad_inches=0.01)

#subprocess.run('mogrify -bordercolor white -trim '+figdir+'*.png',shell=True) 

