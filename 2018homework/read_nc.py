from scipy.io import netcdf
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt

f = netcdf.netcdf_file('~/renql/project/2018homework/ECMWF_ERA-40_subset.nc','r')
        print(f.variables.keys()) #直接获取所有的变量名称
        lat = f.variables['latitude'][:].copy()
        lon = f.variables['longitude'][:].copy()
        tcw = f.variables['tcw'][:].copy()
        print(type(tcw))
        scale_factor=f.variables['tcw'].scale_factor
        add_offset=f.variables['tcw'].add_offset
        print(scale_factor)
        print(add_offset)
        f.close()

        tcw_p=tcw*scale_factor+add_offset
        #print(tcw_p)

        #把X，Y传入网格中
        x,y=np.meshgrid(lon,lat)
        #8:8+2=10,将tcw分为10部分  #alpha：透明度为0.3  
        pic=plt.contourf(x,y,tcw_p[0,:,:],8,cmap=mpl.cm.PuBu,alpha=0.3,linewidth=.5)
        #添加标题：
        plt.title('TCW')
        #坐标轴说明：
        plt.xlabel('Lontituede')
        plt.ylabel('Latitude')
        plt.colorbar(pic)
        plt.show()
