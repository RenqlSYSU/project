#!/bin/sh

#20190309
echo "start copy F2000_CAM5"
scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/project/TP_NUDG/z_tendency-20191022/mdata/NUDG6h-Clim_advetc_dzdt_month.dat /home/ys17-19/renql/project/TP_NUDG/z_tendency-20191022/mdata/ >> /home/ys17-19/renql/project/TP_NUDG/z_tendency-20191022/scp_log
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/F2000_CAM5/post_nudg_5-8/F2000_CAM5.cam.h1* /home/ys17-19/renql/model/F2000_CAM5_NUDG5-8/  >> /home/ys17-19/renql/model/scp_log
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/F2000_CAM5/post_nudg/F2000_CAM5.cam.h1.* /home/ys17-19/renql/model/F2000_CAM5_NUDG/  >> /home/ys17-19/renql/model/scp_log
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/project/TP_NUDG/energy-20180417/1909-draw_clim_month_ave_eddy_qtran_int_shaded_ratio.ncl /home/ys17-19/renql/project/TP_NUDG/energy-20180417/1909-draw_clim_month_ave_eddy_qtran_int_shaded_ratio.ncl >> /home/ys17-19/renql/model/scp_log
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/project/TP_NUDG/analysis/mdata/*2month_ave_intq1.nc /home/ys17-19/renql/project/TP_NUDG/analysis/mdata/  >> /home/ys17-19/renql/model/scp_log

#echo "start copy F2000_F19_CAM4_CTRL_daily"
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/CESM/input/atm/cam/topo/USGS-gtopo30_1.9x2.5_remap_c050602.nc /home/ys17-19/renql/model/TP_CTRL/  >> /home/ys17-19/renql/model/scp_log

#echo "start copy AMIP-CTRL"
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM/post_data/AMIP_C5PM* /home/ys17-19/renql/model/AMIP-CTRL/  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM/post_data/*.V.nc /home/ys17-19/renql/model/AMIP-CTRL  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM/post_data/*.T.nc /home/ys17-19/renql/model/AMIP-CTRL  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM/post_data/*.OMEGA.nc /home/ys17-19/renql/model/AMIP-CTRL >> /home/ys17-19/renql/model/scp_log 

#echo "start copy TP-NUDG-24h"
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_24h/AMIP_C5PM* /home/ys17-19/renql/model/TP-NUDG-24h/  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_24h/*.V.nc /home/ys17-19/renql/model/TP-NUDG-24h  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_24h/*.T.nc /home/ys17-19/renql/model/TP-NUDG-24h  >> /home/ys17-19/renql/model/scp_log 
#scp -r yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_24h/*.OMEGA.nc /home/ys17-19/renql/model/TP-NUDG-24h  >> /home/ys17-19/renql/model/scp_log 

#echo "start copy TP-NUDG-6h"
#scp yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_6h/AMIP_C5PM* /home/ys17-19/renql/model/TP-NUDG-6h/  >> /home/ys17-19/renql/model/scp_log 
#scp yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_6h/*.V.nc /home/ys17-19/renql/model/TP-NUDG-6h/  >> /home/ys17-19/renql/model/scp_log 
#scp yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_6h/*.T.nc /home/ys17-19/renql/model/TP-NUDG-6h/  >> /home/ys17-19/renql/model/scp_log 
#scp yangsong3@hpc4.sysu.edu.cn:/users/yangsong3/renql/F/AMIP_C5PM_TP_NUDG/post_data_6h/*.OMEGA.nc /home/ys17-19/renql/model/TP-NUDG-6h/  >> /home/ys17-19/renql/model/scp_log 

