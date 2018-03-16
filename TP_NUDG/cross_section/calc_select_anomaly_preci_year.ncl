;***************************************************************************
;filename:     calc_select_anomaly_preci_year.ncl 
;read the 27 years PRECC and PRECL of three cases at East Asia(20N£­45N£¬107-130E)
;then calculate the anomalous rainfall var(ncase,nlat,ntime)
;                                            by Ql_Ren
;                                           2018/03/16
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

