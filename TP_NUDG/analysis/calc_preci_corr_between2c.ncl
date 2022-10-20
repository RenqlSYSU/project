;***************************************************************************
;filename:      calc_preci_corr_between2c.ncl
;output_file:   NUDG6h-Preci_corr_with_ctrl.nc 
;read the 27 years PRECC and PRECL of two case at East AsiaI(15S-55N£¬30-150E) 
;and then calculate the correlation between these two cases to get (nmonth,nlat,nlon).
;                                            by Ql_Ren
;                                           2018/01/18
;modification
;2018/02/06 calculate monthly average rainfall first then calculate the correlation
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"NUDG24h","NUDG6h"/) + "-Preci_corr_with_ctrl.nc"
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
;lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var1        = new((/nyear,365,nlat,nlon/),float)
var         = new((/nyear,365,nlat,nlon/),float)
month_var1  = new((/nyear,12 ,nlat,nlon/),float)  ;store the monthly average rainfall
season_var1 = new((/nyear, 4 ,nlat,nlon/),float)  ;store the seasonal average rainfall
month_var   = new((/nyear,12 ,nlat,nlon/),float)  ;store the monthly average rainfall
season_var  = new((/nyear, 4 ,nlat,nlon/),float)  ;store the seasonal average rainfall
month_corr  = new((/      12 ,nlat,nlon/),float)
season_corr = new((/       4 ,nlat,nlon/),float)

var1 = 0 
do nt = 0,nyear-1,1
    do nv = 0,nvar-1,1 
        f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var1(nt,:,:,:) = var1(nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do  ;the loop of PRECC and PRECL
end do      ;the loop of 27 years
var1 = var1*1000*3600*24 ;convert the unit "m/s" to "mm/day"

do nf = 1 , ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
var = 0
do nt = 0,nyear-1,1
    do nv = 0,nvar-1,1 
        f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do  ;the loop of PRECC and PRECL
end do      ;the loop of 27 years
var  = var*1000*3600*24

time  = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    month_var1(:,nm,:,:) = dim_avg_n_Wrap(var1(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
     month_var(:,nm,:,:) =  dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    month_corr(nm,:,:)   = escorc_n(month_var1(:,nm,:,:),month_var(:,nm,:,:),0,0) ;get (nmonth,nlat,nlon)
end do

do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
        select_time := ind(time(:,1).eq.12.or.time(:,1).eq.1.or.time(:,1).eq.2)
    else
        select_time := ind(time(:,1).ge.nm-1.and.time(:,1).le.nm+1)
    end if
    season_var1(:,ns,:,:) = dim_avg_n_Wrap(var1(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    season_var(:,ns,:,:)  = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    season_corr(ns,:,:)   = escorc_n(season_var1(:,ns,:,:),season_var(:,ns,:,:),0,0) ;get (nseason,nlat,nlon)
end do

printVarSummary(month_corr)
printVarSummary(season_corr)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf-1))
ncdf = addfile(fileout(nf-1),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the monthly and seasonal correlation between CTRL precipitation and NUDG"
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
filevardef(ncdf,"month_corr" ,typeof(month_corr) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_corr",typeof(season_corr),(/"season","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_corr   = (/month_corr/)
ncdf->season_corr  = (/season_corr/)
end do
end

