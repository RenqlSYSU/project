;***************************************************************************
;filename:     calc_select_anomaly_preci_year.ncl 
;read the 27 years PRECC and PRECL of three cases at three areas
;of "North China","Central East China" and "South China"
;then calculate the area and JJA average to get the standard deviation and anomalous
;next select the anomaly great than the deviation.
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
  case   = (/"AMIP_CTRL","NUDG24h","NUDG6h"/);this can be used as left string

area  = (/"North_China","Central_East_China","South_China"/)
lats  = (/   34    ,   28    ,   20    /)
latn  = (/   43    ,   34    ,   28    /)
lonl  = (/   107   ,   107   ,   107   /)
lonr  = (/   122   ,   122   ,   122   /)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;==========================================================
f           = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
narea = dimsizes(area)
var   = new((/ncase,narea,nyear/),float)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
var = 0
do nc = 0, ncase-1,1
do nt = 0, nyear-1,1
    do nv = 0, nvar-1,1 
    do na = 0, narea-1,1
        f     = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
        temp := f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        temp := dim_avg_n_Wrap(temp,0)
        var(nc,na,nt) = var(nc,na,nt) + wgt_areaave_Wrap(temp,1.0,1.0,0)
    end do  ;the loop of narea
    end do  ;the loop of PRECC and PRECL
end do   ;the loop of 27 years
end do   ;the loop of 3 cases 
var  = var*1000*3600*24  ;convert unit from m/s to mm/day
var  = dim_standardize_n_Wrap(var,1,2) ;get the standardized anomalies

do nc = 0, ncase-1,1
    do na = 0, narea-1,1
        nega  := ind(var(nc,na,:).lt.-1)
        posi  := ind(var(nc,na,:).gt.1)
        print(case(nc) + " " +area(na) + " negative years: " + years(nega)) 
        print(case(nc) + " " +area(na) + " positive years: " + years(posi))
    end do
end do
end


