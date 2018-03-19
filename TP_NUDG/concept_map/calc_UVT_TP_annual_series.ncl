;***************************************************************************
;filename:      calc_UVT_TP_annual_series.ncl
;                                            by Ql_Ren
;                                           2017/03/19
;******************************************************************************
begin
;the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/concept_map/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-UVT_TP_annual_series.nc"
case     = (/"AMIP_CTRL","NUDG24h","NUDG6h"/)
var_name = (/"U","V","T"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of TP 
lats = 25 
latn = 40
lonl = 70
lonr = 105
lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
var         = new((/nvar,nyear,365,nlev/),float)

do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nv = 0,nvar-1,1
do nt = 0,nyear-1,1
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    temp = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    var(nv,nt,:,:) = wgt_areaave_Wrap(temp,1.0,1.0,0)
end do  ;the loop of the 27 years
end do  ;the loop of the 3 variables
ave = dim_avg_n_Wrap(var,1) ;get the climate (nvar,ndays,nlev)
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
fileAtt@discription   = "The annual cycle of TP average's UVT 850,500,200 hPai,unit is m/s,m/s,K" 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"var","year","day","lev"/)
dimSizes = (/nvar , nyear,365  ,nlev /)
dimUnlim = (/False, False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var" ,"integer",(/"var"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day" ,"integer",(/"day"/))
filevardef(ncdf,"lev" ,"integer",(/"lev"/))
filevardef(ncdf,"series",typeof(var),(/"var","year","day","lev"/))
filevardef(ncdf,"ave"   ,typeof(ave),(/"var","day","lev"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"series",var)
filevarattdef(ncdf,"ave",ave)

ncdf->var    = (/1,2,3/)  ;the three variables are U,V,T
ncdf->year   = (/years/)
ncdf->day    = ispan(1,365,1)
ncdf->lev   = (/lev/)
ncdf->series = (/var/)
ncdf->ave = (/ave/)
end do
end

