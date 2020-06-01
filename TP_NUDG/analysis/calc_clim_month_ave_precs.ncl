;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;read the 27 years variable(TLLL) of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
filein  = (/"/home/ys17-19/data-observation/renql/AMIP_C5PM/pre_data/AMIP_C5PM.cam.h0.",\
            "/home/ys17-19/data-observation/renql/AMIP_C5PM_TP_NUDG/pre_data_6h/AMIP_C5PM_TP_NUDG.cam.h0.",\
            "/home/ys17-19/data-observation/renql/AMIP_C5PM_TP_NUDG/pre_data_24h/AMIP_C5PM_TP_NUDG.cam.h0."/)
fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_ave_precs"
var_name  = (/"PRECSC","PRECSL"/)
years    = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats =  0  ;-20
latn = 90  ;60
lonl =  0  ;20
lonr = 360 ;220

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f1    = addfile(filein(0)+"1995-01.nc","r")
vars:= f1->$var_name(0)$(:,{lats:latn},{lonl:lonr})
vars@units = "mm/day"
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var        = new((/nyear*12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)

do nf = 0 , 1, 1 ;ncase-1,1
    files  = systemfunc("ls " + filein(nf) + "*" + ".nc")
    setfileoption("nc","SuppressClose",False) ; close too many files
    f      = addfiles(files(0:(nyear*12-1)),"r")
    ListSetType(f,"cat")
    time = cd_calendar(f[:]->time,0)
    
    var = 0.0
    do nv = 0, nvar-1,1
    var = var + f[:]->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do
    var  = var*1000*3600*24 ;convert m/s to mm/day

;----------calculate the month and season average-----------------------
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    month_ave(nm,:,:) = dim_avg_n_Wrap(var(select_time,:,:),0) ;get (nyear,nmonth,nlat,nlon)
    month_vari(nm,:,:) = dim_variance_n_Wrap(var(select_time,:,:),0) ;get (nyear,nmonth,nlat,nlon)
end do

do nt = 0,nyear-1,1
do ns = 0,3,1
    nm = (nt*4+ns)*3
    if(ns.eq.0.and.nt.eq.0) then
    var_season(nt,ns,:,:) = (var(0,:,:) + var(1,:,:) + var(nyear*12-1,:,:))/3
    else
    var_season(nt,ns,:,:) = (var(nm-1,:,:) + var(nm,:,:) + var(nm+1,:,:))/3
    end if
end do
end do
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
season_vari = dim_variance_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
printVarSummary(month_ave)
printVarSummary(season_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+".nc")
ncdf = addfile(fileout(nf)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the monthly and seasonal average"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","season","lat","lon"/)
dimSizes = (/  12   ,  4     ,nlat ,nlon/)
dimUnlim = (/False  ,False   ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",vars)
filevarattdef(ncdf,"season_ave",vars)
filevarattdef(ncdf,"month_vari",vars)
filevarattdef(ncdf,"season_vari",vars)

ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end


