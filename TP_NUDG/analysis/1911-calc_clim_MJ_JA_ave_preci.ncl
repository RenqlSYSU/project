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
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_MJ_JA_ave_preci.nc"
    years    = ispan(1979,2005,1)   ;the calculate time
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_preci.nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

var_name = (/"PRECC","PRECL"/) ;,"OMEGA","RELHUM"

;area: east Asia
lats = -20 ;-90 ;
latn = 60  ; 90 ;
lonl = 20  ; 0  ;
lonr = 220 ;360 ;
month_s = 5
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
time1 = cd_calendar(f->time,0)
select_time1 = ind(time1(:,1).ge.month_s.and.time1(:,1).le.month_e)
nmonth = 2
ntime = dimsizes(select_time1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
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
var = 0
do nv = 0,nvar-1,1 
do nt = 0,nyear-1,1
    f   = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(select_time1,{lats:latn},{lonl:lonr})
end do
end do
var  = var*1000*3600*24 ;convert m/s to mm/day

time = cd_calendar(f->time(select_time1),0)
select_time := ind(time(:,1).ge.5.and.time(:,1).le.6)
var_month(:,0,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
select_time := ind(time(:,1).ge.7.and.time(:,1).le.8)
var_month(:,1,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)

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
fileAtt@discription   = "ealy(56) and late(78) stage of the summer"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","lat","lon"/)
dimSizes = (/  2    ,nlat ,nlon/)
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

ncdf->month  = (/56,78/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end do
;end do ;loop of variable
end

