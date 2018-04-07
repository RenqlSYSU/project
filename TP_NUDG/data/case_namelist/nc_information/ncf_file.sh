#!/bin/sh
source ~/.bashrc
#ncf /users/yangsong3/CESM/input/atm/cam/solar/SOLAR_SPECTRAL_Lean_1610-2008_annual_c090324.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/topo/USGS-gtopo30_0.9x1.25_remap_c051027.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/physprops/water_refindex_rrtmg_c080910.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/physprops/iceoptics_c080917.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/physprops/F_nwvl200_mu20_lam50_res64_t298_c080428.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/volc/CCSM4_volcanic_1850-2008_prototype1.nc >> atm_in_ncfile
#ncf /users/yangsong3/CESM/input/atm/cam/chem/trop_mozart/ub/clim_p_trop.nc >> atm_in_ncfile
#ncf  >> atm_in_ncfile

ncf /users/yangsong3/CESM/input/atm/cam/sst/sst_HadOIBl_bc_0.9x1.25_1850_2012_c130411.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/pftdata/pft-physiology.clm40.c130424.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/surfdata/surfdata.pftdyn_0.9x1.25_simyr1850-2005_c091008.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/snicardata/snicar_optics_5bnd_c090915.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/surfdata/surfdata_0.9x1.25_simyr1850_c110921.nc >> docn_clm_in_ncfile
ncf /users/yangsong3/CESM/input/lnd/clm2/rtmdata/rdirc_0.5x0.5_simyr2000_slpmxvl_c120717.nc >> docn_clm_in_ncfile

ncf /users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/b40_20th_1d_b08c5cn_139jp.clm2.r.1979-01-01-00000.nc >> initial_ncfile
ncf /users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/b40_20th_1d_b08c5cn_139jp.cice.r.1979-01-01-00000.nc >> initial_ncfile
ncf /users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/b40_20th_1d_b08c5cn_139jp.cam.i.1979-01-01-00000.nc >> initial_ncfile
#ncf  >> initial_ncfile

