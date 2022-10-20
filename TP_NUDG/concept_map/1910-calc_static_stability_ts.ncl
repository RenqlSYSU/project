;***************************************************************************
;filename:      calc_UVT_TP_annual_series.ncl
;                                            by Ql_Ren
;                                           2017/03/19
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case    = (/"CTRL","NUDG6h","NUDG24h"/);"NG24h",
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/concept_map/mdata/" + case + "-static_2EA_annual_series.nc" 
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_EKE-" + cb + "-" + ca + "day.nc"
end if

var_name = (/"T","Q"/)

area = (/"SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /); 25  ,25   ,30   , ;-20
latn = (/ 30  , 45  /); 40  ,30   ,35   , ;60
lonl = (/107  ,115  /); 70  ,70   ,70   , ;20
lonr = (/130  ,130  /); 105 ,105  ,105  , ;130 ;220
lev  = (/925,600/)
narea = dimsizes(area)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(case)
var   = new((/narea,nyear,365/),float)

do nc = 0, ncase-1,1
do na = 0,narea-1,1
print(area(na))
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
q     := new((/nyear,ntime,nlev,nlat,nlon/),float)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nt = 0, nyear-1 ,1
    nv = 0
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    t(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nv = 1
    f    = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    q(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do  ;the loop of the 27 years
theta := pot_temp_equiv(lev*100,t,q,2,"q") ;Compute the equivalent potential temperature
var(na,:,:) = wgt_areaave((theta(:,:,0,:,:)-theta(:,:,1,:,:))/3.25,1.0,1.0,0)
end do

var@units = "K/100hPa"
var@long_name = " the vertical gradient of equivalent potential temperature between 925 and 600 hPa"
ave = dim_avg_n_Wrap(var,1) 
printVarSummary(theta)
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
fileAtt@discription   = "The annual cycle of SEA(" + lats(0) + "-" + latn(0) + "N," + lonl(0) + "-" + lonr(0) + \
                        "E) and NEA(" + lats(1) + "-" + latn(1) + "N," + lonl(1) + "-" + lonr(1) + "E)"
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


