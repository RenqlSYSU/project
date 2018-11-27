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
fileout = path + "/concept_map/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-preci_6TP_annual_series.nc"
case     = (/"CTRL","NG24h","NG6h"/)
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

area = (/"SETP","SCTP","SWTP","CETP","CWTP","NTP"/);"ESTP","ECTP","NTP","SEA","NEA"/) "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/25    ,25    ,25    ,30    ,30    ,35   /);25    ,30    ,35   , 22  , 30  /) 25  ,25   ,30   , ;-20
latn = (/30    ,30    ,30    ,35    ,35    ,40   /);30    ,35    ,40   , 30  , 45  /) 40  ,30   ,35   , ;60
lonl = (/100   ,85    ,70    ,90    ,70    ,70   /);100   ,90    ,70   ,107  ,107  /) 70  ,70   ,70   , ;20
lonr = (/105   ,100   ,85    ,105   ,90    ,105  /);105   ,105   ,105  ,130  ,130  /) 105 ,105  ,105  , ;130 ;220

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
narea = dimsizes(area)
var   = new((/narea,nyear,365/),float)

do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
var = 0
do na = 0,narea-1,1
do nt = 0,nyear-1,1
do nv = 0,nvar-1,1
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    temp := f->$var_name(nv)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
    var(na,nt,:) = var(na,nt,:) + wgt_areaave_Wrap(temp,1.0,1.0,0)
end do  ;the loop of the 27 years
end do  ;the loop of the 3 variables
end do  
var = var*1000*24*3600
var@units = "mm/day"
ave = dim_avg_n_Wrap(var,1) ;get the climate (narea,ndays)
printVarSummary(var)
printVarSummary(ave)
printVarSummary(temp&time)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The annual cycle of area average's preci,unit is mm/day" + \
        "the area is SETP,SCTP,SWTP,CETP,CWTP,NTP"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"area","year","day"/)
dimSizes = (/ narea,nyear,365  /)
dimUnlim = (/ False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"area","integer",(/"area"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day" ,"integer",(/"day"/))
filevardef(ncdf,"series",typeof(var),(/"area","year","day"/))
filevardef(ncdf,"ave"   ,typeof(ave),(/"area","day"/))

;Define the attribute of the variables
filevarattdef(ncdf,"series",var)
filevarattdef(ncdf,"ave",ave)

ncdf->area   = ispan(1,narea,1)
ncdf->year   = (/years/)
ncdf->day    = ispan(1,365,1)
ncdf->series = (/var/)
ncdf->ave = (/ave/)
end do
end

