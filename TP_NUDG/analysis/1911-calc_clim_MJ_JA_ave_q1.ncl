;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;read the 27 years variable(TLLL) of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
;
;Modification:
;20190712: add the case of noTP
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_MJ_JA_ave_Q1.nc"
    years    = ispan(1979,2005,1)   ;the calculate time
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_q1.nc"
    years   = ispan(21,30,1)   ;the calculate time
end if
if(case_option.eq.2) then
    case    = (/"F2000","NUDG"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_Q1.nc"
end if

var_name= (/"Q1","PS"/)
g  = 9.8 ;m/(s*s)

;area: east Asia
lats =  0  ;-90 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;
lev  = (/1000,925,850,700,600,500,400,300,250,200,150,100/) ;hPa
;lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
month_s = 1
month_e = 12
nmonth = 12

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time1 = cd_calendar(f->time,0)
select_time1 = ind(time1(:,1).ge.month_s.and.time1(:,1).le.month_e)
ntime = dimsizes(select_time1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var        = new((/nyear,ntime ,nlat,nlon/),float)
var_month  = new((/nyear,nmonth,nlat,nlon/),float)
month_ave  = new((/      nmonth,nlat,nlon/),float)
month_vari = new((/      nmonth,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups

do nf = 0 , ncase-1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
q1    = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("handle with " + years(nt))
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q1(nt,:,:,:,:)  = tofloat(f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr}))

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
ps(nt,:,:,:)   = f->$var_name(nv)$(select_time1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
print(q1(10,10,2,10:20,30))
var = dim_sum_n(q1*dp,2)/g
delete([/q1,ps,dp/])

var@long_name = "Integrated Apparent Heat Source calc from temperature equation"
var@units     = "W/kg" 

time = cd_calendar(f->time(select_time1),0)
do nm = 0, nmonth-1, 1
mons = nm + 1
mone = nm + 1
select_time := ind(time(:,1).ge.mons.and.time(:,1).le.mone)
var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
copy_VarMeta(var(0,0,:,:),var_month(0,0,:,:))
printVarSummary(month_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "Integrated Apparent Heat Source calculated from temperature equation"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","lat","lon"/)
dimSizes = (/nmonth ,nlat ,nlon/)
dimUnlim = (/False  ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"month_vari",typeof(month_vari),(/"month" ,"lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var_month)
filevarattdef(ncdf,"month_vari",var_month)

ncdf->month  = ispan(1,12,1) ;(/56,78/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end do
;end do ;loop of variable
end

