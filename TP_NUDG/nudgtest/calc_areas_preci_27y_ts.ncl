;***************************************************************************
;filename:      calc_areas_preci_27y_ts.ncl
;output_file:   Areas_preci_27y_ts.nc 
;read the 27 years preci(PRECC+PRECL) in one case and then 
;calculate the monthly and seasonly time series of seven area in order to test the right of the 
;NUDG experiment.
;                                            by Ql_Ren
;                                           2017/12/21
;******************************************************************************
begin
;the path of the file
;filein   = "/users/yangsong3/renql/project/TP_NUDG/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.";+year+".daily.nc"
;filein   = "/users/yangsong3/renql/project/TP_NUDG/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR." ;+year+".daily.nc"
;fileout  = "/users/yangsong3/renql/project/TP_NUDG/nudgtest/mdata/CTRL-Areas_preci_27y_ts.nc"
;fileout  = "/users/yangsong3/renql/project/TP_NUDG/nudgtest/mdata/NUDG24h-Areas_preci_27y_ts.nc"
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/nudgtest/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Areas_preci_27y_ts.nc"

var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"India","Tibet","Philippines","Maritime_continent",\
          "North China","Central East China","South China","TPnudg"/)
lats = (/20 ,25  ,5   ,-10 ,34  ,28  ,20 ,25 /)
latn = (/35 ,35  ,20  ,7.5 ,43  ,34  ,28 ,40 /)
lonl = (/60 ,90  ,120 ,85  ,107 ,107 ,107,70 /)
lonr = (/80 ,103 ,140 ,120 ,122 ,122 ,122,105/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
nyear = dimsizes(years)
narea = dimsizes(area)
ncase = dimsizes(filein)
index = new((/narea,nyear,365/),float)
month_series  = new((/narea,nyear,12/),float)
season_series = new((/narea,nyear,4/),float)

do nf = 0, ncase-1,1
;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
index = 0
do na = 0,narea-1,1
    do nt = 0,nyear-1,1
        do v = 0, dimsizes(var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
            f  = addfile(filein(nf)+years(nt)+".daily."+var_name(v)+".nc","r")
            vars := f->$var_name(v)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
            index(na,nt,:) = index(na,nt,:) + wgt_areaave_Wrap(vars,1.0,1.0,0)*3600*24*1000 ;calculate the area average and get the index
        end do ;the loop of precc and precl
    end do  ;the loop of the years
end do   ;the loop of the areas
printVarSummary(index)

;---------calculate the month average------------------------
time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    month_series(:,:,nm) = dim_avg_n_Wrap(index(:,:,select_time),2)  ;get var(narea,nyear,12)
end do

;----------calculate the season average-----------------------
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    season_series(:,:,ns) = ( month_series(:,:,11)*mdays(11) + \
                    month_series(:,:,nm)*mdays(nm) + month_series(:,:,nm+1)*mdays(nm+1))\
                    /(mdays(11) + mdays(nm)+ mdays(nm+1))
    else
    season_series(:,:,ns) = ( month_series(:,:,nm-1)*mdays(nm-1) + \
                    month_series(:,:,nm)*mdays(nm) + month_series(:,:,nm+1)*mdays(nm+1))\
                    /(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
    end if
end do
;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "monthly and seasonly preci 27years time series of seven areas " + \
                    "the seven areas are India,Tibet,Philippines,Maritime_continent," +\
                    "North China,Central East China,South China"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"area","year","month","season"/)
dimSizes = (/narea,nyear,12,4/)
dimUnlim = (/False, False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"area","integer",(/"area"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"month_series",typeof(month_series),(/"area","year","month"/))
filevardef(ncdf,"season_series",typeof(season_series),(/"area","year","season"/))

;Define the attribute of the variables
month_series@long_name = "the 27 years time series of seven areas preci"
month_series@units     = "mm/day"
season_series@long_name = "the 27 years time series of seven areas preci"
season_series@units     = "mm/day"
filevarattdef(ncdf,"month_series",month_series)
filevarattdef(ncdf,"season_series",season_series)

ncdf->area   = (/1,2,3,4,5,6,7,8/)  ;the two cases are (/"CTRL","NUDG"/)
ncdf->year   = (/years/)
ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->month_series = (/month_series/)
ncdf->season_series = (/season_series/)
end do
end


