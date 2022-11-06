#!/usr/bin/env python
import cf
import cfplot as cfp
import numpy as np
import sys
import subprocess
import matplotlib

f0=cf.read("/home/users/qd201969/gtopo30_0.9x1.25.nc")
phis=f0[2]
print(repr(phis))
phis=phis/9.8 # transfer from m2/s2 to m

path = '/home/users/qd201969/ERA5-1HR-lev/'
#filename = ['ff_trs_pos','tr_trs_pos']
lev = ['850','500','250']
nl1 = 2
nl2 = 1
filename = [lev[nl1]+'_'+lev[nl2]+'_no',lev[nl1]+'_'+lev[nl2]+'_yes' \
           ,lev[nl2]+'_'+lev[nl1]+'_no',lev[nl2]+'_'+lev[nl1]+'_yes']
#case='match'+lev[nl1]+'_'+lev[nl2]
case='pas_match'+lev[nl1]+'_'+lev[nl2]
filt = False #True #

#f=cf.read('/home/users/qd201969/TRACK_result/ERA5_VOR850_6hr_2000_DET_T42/ff_trs_pos.nc')
#f=cf.read('/home/users/qd201969/ERA5-1HR/ff_trs_pos.oct-mar1979-2020_addmslp_addwind925_addwind10m.new.nc')
#f=cf.read('/home/users/qd201969/ERA5-1HR/ff_trs_pos.apr-sep1979-2019_addmslp_addwind925_addwind10m.new.nc')

for nc in range(0,len(filename),1):#,len(f),1):
    figname  = str(case+'_'+filename[nc]+'.png')
    #f=cf.read(path+'ERA5_VOR'+str(lev[nl])+'_6hr_2000_DET_T42/'+filename[nc]+'.nc')
    f=cf.read(path+case+'/'+filename[nc]+'.nc')
    print(f)
    #g = f[11]
    g = f[2]
    g = g*1e5
    print(g)

    if filt:
        x = g.aux('longitude').data
        y = g.aux('latitude').data
        tim = g.aux('time').data

        nt = 0
        ni = range(0,g.shape[0])
        nj = range(0,g.shape[1])

        for i in ni:
            inreg = 0
            h = g[i].where(g[i].array>1e+10, cf.masked)
            print(h.shape)
            h = np.squeeze(h.array)
            print(h.shape)

            for j in nj:
                if x[i,j] >= 30 and x[i,j] <= 110 and y[i,j] >= 20 and y[i,j] <= 55 and not x[i,j].mask and not y[i,j].mask:
                    jmax = np.argmax(h)
                    tmax = h[jmax]
                    longmax = x[i, jmax]
                    latmax  = y[i, jmax]
                    print(longmax, latmax, tmax)
                    inreg = 1
                    break

            if inreg == 1:
                if nt == 0:
                    tr_arb = g[i]
                    nt = 1
                else:
                    tr_arb = [tr_arb, g[i]]

        tr_arb = cf.aggregate(tr_arb, relaxed_identities=True)
        g = tr_arb[0]

    cfp.setvars(file=figname)
    cfp.gopen()
    cfp.mapset(lonmin=0, lonmax=150, latmin=10, latmax=70)

    cfp.levs(manual=[1500,3000,4500])
    cfp.con(phis,fill=False, lines=True, colors='k',linewidths=2)
    cfp.levs()

    cfp.cscale('precip2_17lev')
    cfp.levs(min=0.0, max=8.0, step=0.5)
    #cfp.traj(g, linewidth=2, legend_lines=False, markeredgecolor='none', markeredgewidth=0.0, markevery=1, legend=True, markersize=2)
    #cfp.traj(g)
    cfp.traj(g, zorder=0, legend_lines=True,  linewidth=2, title=filename[nc],colorbar_title=' Relative Vorticity (Hz) * 1e5')
    
    cfp.gclose()

subprocess.run('mogrify -bordercolor white -trim ./'+case+'*.png',shell=True) 

