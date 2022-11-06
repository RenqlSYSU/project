import cf
import cfplot as cfp
import subprocess
import matplotlib
matplotlib.use('Agg')

lonl=0
lonr=360
lats=0
latn=90

lev=[[0    ,1600 ,100 ], # 0Feature Density
     [0    ,32   ,2   ], # 1Genesis Density
     [0    ,32   ,2   ], # 2Lysis Density
     [-1   ,1    ,1   ], # 3Mean Area
     [-1.6 ,1.6  ,0.2 ], # 4Mean Growth/Decay Rate
     [-1   ,1    ,1   ], # 5Mean Anisotropy
     [0    ,8    ,0.5 ], # 6Mean Lifetime
     [0    ,80   ,10  ], # 7Mean Speed
     [0    ,8    ,0.5 ], # 8Mean Intensity
     [-0.04,0.04 ,0.005], # 9Mean Tendency
     [-1   ,1    ,1   ], # 10Spare1
     [-1   ,1    ,1   ], # 11Spare2
     [0    ,160  ,10  ], # 12Std of Speed
     [0    ,6.4  ,0.4 ], # 13Std of Intensity
     [0    ,160  ,10  ], # 14Track Density
     [-1   ,1    ,1   ], # 15X-component of Mean Orientation Vector
     [-80  ,80   ,10  ], # 16X-component of Mean Velocity
     [-1   ,1    ,1   ], # 17Y-component of Mean Orientation Vector
     [-80  ,80   ,10  ]] # 18Y-component of Mean Velocity

draw=[0,1,2,14,7,8]
#draw=[4,6,16,18,12,13]

f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis=f0[2]
print(repr(phis))
phis=phis/9.8 # transfer from m2/s2 to m

files=['/home/users/qd201969/ERA5-1HR/stat_oct-mar1979-2020.nc',
       '/home/users/qd201969/ERA5-1HR/stat_apr-sep1979-2019.nc',
       '/home/users/qd201969/ERA5-1HR/stat_clim_1.nc']
titls=['ONDJFM','AMJJAS','Jan']
nc = 2
    
f=cf.read(files[nc])

cfp.setvars(file='stat_'+titls[nc]+'.png')
cfp.gopen(figsize=[20, 20],rows=3,columns=2,wspace=0.1,hspace=0.001,bottom=0.45)
cfp.mapset(lonmin=lonl, lonmax=lonr, latmin=lats, latmax=latn)

for nv in range(0,len(draw),1):#,len(f),1):
    np = nv+1
    cfp.gpos(np)
    if draw[nv] == 10 or draw[nv] == 11:
        continue

    g=f[draw[nv]]
    if draw[nv] == 4:
        g=g*24

    #g=g.subspace(longitude=cf.wi(lonl,lonr),latitude=cf.wi(lats,latn))
    print(repr(g))

    cfp.levs(min=lev[draw[nv]][0], max=lev[draw[nv]][1], step=lev[draw[nv]][2])
    if lev[draw[nv]][0] < 0 :
        cfp.cscale('BlueDarkRed18')
    else:
        cfp.cscale('precip2_17lev')

    cfp.con(g,fill=True,lines=False,colorbar_title=titls[nc]+' '+g.long_name)
    cfp.levs(manual=[1500,3000])
    cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
    cfp.levs()
    
cfp.gclose()
subprocess.run('mogrify -bordercolor white -trim /home/users/qd201969/ERA5-1HR/stat_'+titls[nc]+'.png',shell=True) 

