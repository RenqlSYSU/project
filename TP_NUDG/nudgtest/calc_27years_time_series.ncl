;***************************************************************************
;filename:      calc_27years_time_series.ncl
;output_file:   27years_time_series.nc 
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h) and then 
;calculate the season and area average in order to test the right of the 
;NUDG experiment.
;                                            by Ql_Ren
;                                           2017/01/08
;******************************************************************************
begin
;the path of the file
filein   = (/"/users/yangsong3/renql/project/TP_NUDG/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/renql/project/TP_NUDG/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/renql/project/TP_NUDG/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/renql/project/TP_NUDG/nudgtest/mdata/27years_time_series.nc"
var_name = (/"U","V","T"/)
case     = (/"AMIP_CTRL","TP_NUDG-6h","TP_NUDG-24h"/)
years    = ispan(1979,2005,1)   ;the calculate time
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)

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
var   = new((/nvar,ncase,nyear,365,nlev,nlat,nlon/),float)
month_series = new((/nvar,ncase,nyear,12,nlev/),float)
season_series = new((/nvar,ncase,nyear,4,nlev/),float)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
do nv = 0,nvar-1,1 
do nc = 0,ncase-1,1
    do nt = 0,nyear-1,1
        f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nv,nc,nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    end do  ;the loop of the 27 years
end do   ;the loop of the three cases 
end do   ;the loop of the 3 variables
var := wgt_areaave_Wrap(var,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)

;---------calculate the month average------------------------
time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    month_series(:,:,:,nm,:) = dim_avg_n_Wrap(var(:,:,:,select_time,:),3)  ;get var(nvar,ncase,nyear,12,nlev)
end do

;----------calculate the season average-----------------------
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    season_series(:,:,:,ns,:) = ( month_series(:,:,:,11,:)*mdays(11) + \
                    month_series(:,:,:,nm,:)*mdays(nm) + month_series(:,:,:,nm+1,:)*mdays(nm+1))\
                    /(mdays(11) + mdays(nm)+ mdays(nm+1))
    else
    season_series(:,:,:,ns,:) = ( month_series(:,:,:,nm-1,:)*mdays(nm-1) + \
                    month_series(:,:,:,nm,:)*mdays(nm) + month_series(:,:,:,nm+1,:)*mdays(nm+1))\
                    /(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
    end if
end do

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The two case are CTRL,NUDG6h,NUDG24h" + \
                    "and the three variables are U,V,T"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"var","case","year","month","season","lev"/)
dimSizes = (/nvar,ncase,nyear,12,4,nlev/)
dimUnlim = (/False, False,False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var","float",(/"var"/))
filevardef(ncdf,"case","float",(/"case"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"month_series",typeof(month_series),(/"var","case","year","month","lev"/))
filevardef(ncdf,"season_series",typeof(season_series),(/"var","case","year","season","lev"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
month_series@long_name = "the time series of three variables and three case of 27 years"
month_series@units     = ""
season_series@long_name = "the time series of three variables and three case of 27 years"
season_series@units     = ""
filevarattdef(ncdf,"month_series",month_series)
filevarattdef(ncdf,"season_series",season_series)

ncdf->var    = (/1,2,3/)  ;the three variables are U,V,T
ncdf->case   = (/1,2,3/)  ;the two cases are (/"CTRL","NUDG6h","NUDG24h"/)
ncdf->year   = (/years/)
ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->lev   = (/lev/)
ncdf->month_series = (/month_series/)
ncdf->season_series = (/season_series/)
end

