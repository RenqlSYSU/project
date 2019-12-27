;***************************************************************************
;filename:      calc_clim_month_ave_preci.ncl
;output_file:   Clim_month_ave_preci.nc 
;read the 27 years PRECC and PRECL of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
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
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_ave_preci.nc"
    years    = ispan(1979,2005,1)   ;the calculate time
end if

if(case_option.eq.1) then 
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_preci.nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

if(case_option.eq.2) then
    case    = (/"F2000","NUDG","NG58"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_preci.nc"
    years    = ispan(1,30,1)   ;the calculate time
end if

;var_name = (/"SHFLX","LHFLX"/)
var_name = (/"PRECC","PRECL"/)
choice   = 0 ;1 means to calc T change due to Heat Flux ; 0 means to calc precipitation

;area: east Asia
lats = -20 ;-90 ;
latn = 60  ; 90 ;
lonl = 20  ; 0  ;
lonr = 220 ;360 ;
;lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var        = new((/nyear,365,nlat,nlon/),float)
ps         = new((/nyear,365,nlat,nlon/),float)
var_month  = new((/nyear,12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float)
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)
printVarSummary(vars)

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
var = 0
do nv = 0,nvar-1,1 
files = systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
do nt = 0,nyear-1,1
    if(case_option.eq.2) then 
    f   = addfile(files(nt),"r")
    else
    f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    end if
    var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years 
end do  ;the loop of PRECC and PRECL

if(choice.eq.1) then 
;    do nt = 0,nyear-1,1
;        f  = addfile(filein(nf)+years(nt)+".daily.PS.nc","r")
;        ps(nt,:,:,:) = ps(nt,:,:,:) + f->PS(:,{lats:latn},{lonl:lonr})  ;surface pressure
;    end do   ;the loop of 27 years 
    cp = 1004.5 ;specific heat capacities of air, unit is J/(kg*K)
    g  = 9.8 ;gravitational acceleration, N/kg, PS/g is th weight of air, kg
    dst = 0.39 ;air density, unit is kg/m3
    var = var/(cp*dst)
else
    var  = var*1000*3600*24 ;convert m/s to mm/day
end if

time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

;----------calculate the season average-----------------------
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    var_season(:,ns,:,:) = ( var_month(:,11  ,:,:)*mdays(11)   + \
                             var_month(:,nm  ,:,:)*mdays(nm)   + \
                             var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(11) + mdays(nm)+ mdays(nm+1))
    else
    var_season(:,ns,:,:) = ( var_month(:,nm-1,:,:)*mdays(nm-1) + \
                             var_month(:,nm  ,:,:)*mdays(nm)   + \
                             var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
    end if
end do
month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
season_vari = dim_variance_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
printVarSummary(month_ave)
printVarSummary(season_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
;fileAtt@discription   = "the monthly and seasonal average Asia Heat flux climatology"
fileAtt@discription   = "the monthly and seasonal average Temperature change due to Heat flux"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","season","lat","lon"/)
dimSizes = (/  12   ,  4     ,nlat ,nlon/)
dimUnlim = (/False  ,False   ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave",typeof(season_ave),(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
;month_ave@long_name = "precipitation rate"
;month_ave@units     = "mm/day"
filevarattdef(ncdf,"month_ave",month_ave)
;season_ave@long_name = "precipitation rate"
;season_ave@units     = "mm/day"
filevarattdef(ncdf,"season_ave",season_ave)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end

