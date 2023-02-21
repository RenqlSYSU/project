#!/bin/sh
cd /WORK/sysu_hjkx_ys/renql/NUDG
yhrun -n 1 -N 1 -c 1 -p work -J extract sh ./lee-cam_pkg_continue_run_daily.sh > inform.txt

#cd /HOME/sysu_hjkx_ys/WORKSPACE/renql/F/data/
# through yhbatch.sh to submit this file
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./sf_djf.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p6x2_clim_perci_ctrl.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p6x2_clim_wind_ctrl.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_preci_sf_dif_DJF.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_preci_sf_ano_DJF.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_qano.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_vano.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_vano_revise.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_qano_revise.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano_decompose.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano_decompose_v.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano_decompose_q.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano_decompose_vqbar.ncl
#yhrun -n 1 -N 1 -c 1 -p work -J index ncl ./draw_1p5x2_clim_qv_ano_decompose_qvbar.ncl
