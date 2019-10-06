;***************************************************************************
;filename:      calc_UVT_TP_annual_series.ncl
;                                            by Ql_Ren
;                                           2017/03/19
;******************************************************************************
begin
;the path of the file
path    = "/home/ys17-19/renql"
filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/project/TP_NUDG/concept_map/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-preci_2EA_annual_series.nc"
case     = (/"CTRL","NG6h","NG24h"/)
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

area = (/"SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /); 25  ,25   ,30   , ;-20
latn = (/ 30  , 45  /); 40  ,30   ,35   , ;60
lonl = (/107  ,115  /); 70  ,70   ,70   , ;20
lonr = (/130  ,130  /); 105 ,105  ,105  , ;130 ;220
lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
narea = dimsizes(area)
var   = new((/narea,nyear,365,nlev/),float)

do nv = 4,nvar-1,1
fileout = path + "/concept_map/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-" + var_name(nv) + "_TP_annual_series.nc"
do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    temp = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    var(nt,:,:) = wgt_areaave_Wrap(temp,1.0,1.0,0)
end do  ;the loop of the 27 years
ave = dim_avg_n_Wrap(var,0) ;get the climate (nvar,ndays,nlev)
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
fileAtt@discription   = "The annual cycle of TP average" 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"year","day","lev"/)
dimSizes = (/ nyear,365  ,nlev /)
dimUnlim = (/ False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day" ,"integer",(/"day"/))
filevardef(ncdf,"lev" ,"integer",(/"lev"/))
filevardef(ncdf,"series",typeof(var),(/"year","day","lev"/))
filevardef(ncdf,"ave"   ,typeof(ave),(/"day","lev"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",var&level)
filevarattdef(ncdf,"series",var)
filevarattdef(ncdf,"ave",ave)

ncdf->year   = (/years/)
ncdf->day    = ispan(1,365,1)
ncdf->lev   = (/lev/)
ncdf->series = (/var/)
ncdf->ave = (/ave/)
end do
end do
end

