#!/usr/bin/env python
import matplotlib 
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

plt.rcParams["font.weight"] = "bold"
#plt.rcParams['axes.linewidth'] = 3
font = {'family': 'serif',
        'style': 'normal',
        'weight': 'bold',
        'color':  'black', 
        'size': 14,
        }

path = "/home/ys17-19/renql/project/2020MiddleEastJet/"
area = ["ME1","ME2","EA1","EA2","NA1","NA2"]
forc = ["Qd","lQte","lFte","hQte","hFte","Tadv","vort"]
corr = np.loadtxt(path+"project_dzdt_corr.txt")
prob = np.loadtxt(path+"project_dzdt_prob.txt")
print(corr)
print(prob)

def draw_relation(corr,prob,labelx,labely,figdir):
    midfont=10
    smfont=14
    cnlevels = [-0.99,-0.95,-0.90,0.90,0.95,0.99]
    listcolor = ["dodgerblue","deepskyblue","powderblue","white","gold","darkorange","red"]
    fcolors = colors.ListedColormap(listcolor)
    norm = colors.BoundaryNorm(boundaries=cnlevels, ncolors=fcolors.N,extend='both')
    
    fig = plt.figure(figsize=(9,9),dpi=200)
    axe = fig.add_axes([0.05, 0.05, 0.9, 0.6])
    #axe.set_title(title,fontsize=midfont) 
            
    im = axe.imshow(prob, cmap=fcolors, norm=norm)
    axe.set_xticks(np.arange(len(labelx)))
    axe.set_yticks(np.arange(len(labely)))
    axe.set_xticklabels(labelx,fontdict=font)#,fontsize=smfont)
    axe.set_yticklabels(labely,fontdict=font)#,fontsize=smfont)
    axe.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    axe.plot([1.5,1.5],[-0.5,(len(labely)-0.5)],color="k",linewidth=2)
    axe.plot([3.5,3.5],[-0.5,(len(labely)-0.5)],color="k",linewidth=2)

    #axe.spines[:].set_visible(False)
    #axe.spines[0:3].set_linewidth('3.0')
    # Loop over data dimensions and create text annotations.
    for i in range(len(labely)):
        for j in range(len(labelx)):
            text = axe.text(j,i,"%.2f"%corr[i,j], ha="center", va="center",fontdict=font)#, color="k",fontsize=smfont)

    #position = fig.add_axes([0.92, 0.2, 0.01, 0.6]) #left, bottom, width, height
    #cb = plt.colorbar(im, cax=position ,orientation='vertical')#, shrink=.9)
    cb = plt.colorbar(im, ax=axe, shrink=.9,extendfrac='auto')
    cb.ax.tick_params(labelsize=smfont)
    #fig.tight_layout(rect=(0,bmlo,1,1)) #w_pad=0.5,h_pad=0.001) #,
    #fig.savefig(figdir)
    fig.savefig(figdir+".png", bbox_inches='tight',pad_inches=0.08)

prob = 1-prob
prob = np.where(corr<0,-prob,prob)
draw_relation(corr,prob,area,forc,path+"fig/eof_corr")

