;***************************************************************************
;filename:      calc_preci_TP_annual_series.ncl
;
;1 take the sum of PRECC plus PRECL, then convert unit from m/s to mm/day
;2 get the var(nyear,nday) and the climate mean state (nday)
;3 save the data in three files
;                                            by Ql_Ren
;                                           2017/03/23
;******************************************************************************
begin
;the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/concept_map/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-preci_TP_annual_series.nc"
case     = (/"AMIP_CTRL","NUDG24h","NUDG6h"/)
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
var         = new((/nyear,365/),float)

do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
var = 0
do nt = 0,nyear-1,1
do nv = 0,nvar-1,1
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    temp = f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    var(nt,:) = var(nt,:) + wgt_areaave_Wrap(temp,1.0,1.0,0)
end do  ;the loop of the 27 years
end do  ;the loop of the 3 variables
var = var*1000*24*3600
ave = dim_avg_n_Wrap(var,0) ;get the climate (ndays
printVarSummary(var)
printVarSummary(ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The annual cycle of TP average's preci,unit is mm/day" 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"year","day"/)
dimSizes = (/ nyear,365  /)
dimUnlim = (/ False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day" ,"integer",(/"day"/))
filevardef(ncdf,"series",typeof(var),(/"year","day"/))
filevardef(ncdf,"ave"   ,typeof(ave),(/"day"/))

;Define the attribute of the variables
filevarattdef(ncdf,"series",var)
filevarattdef(ncdf,"ave",ave)

ncdf->year   = (/years/)
ncdf->day    = ispan(1,365,1)
ncdf->series = (/var/)
ncdf->ave = (/ave/)
end do
end

