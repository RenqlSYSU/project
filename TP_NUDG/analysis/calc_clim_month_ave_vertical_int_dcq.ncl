;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;read the 27 years variable(TLLL) of one case at East AsiaI(15S��55N��30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = (/"/users/yangsong3/L_Zealot/F/AMIP_C5PM/exe/AMIP_C5PM.cam.h0.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/pre_data_24h/AMIP_C5PM_TP_NUDG.cam.h0.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/AMIP_C5PM_TP_NUDG.cam.h0."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_ave_vertical_int_"
var_name  = (/"DCQ","DTCOND"/)
years    = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f1    = addfile(filein(0)+"1995-01.nc","r")
vars:= f1->$var_name(0)$(:,:,{lats:latn},{lonl:lonr})
hyai = f1->hyai  ;hyai(ilev),ilev = lev+1,hybrid A coefficient at layer interfaces
hybi = f1->hybi
p0   = 100000   ;a scalar ,double,reference pressure,unit is Pa
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(vars&lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var_temp   = new((/nyear*12 ,nlev,nlat,nlon/),float)
var        = new((/nyear*12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)

do nf = 0 , ncase-1,1
    files  = systemfunc("ls " + filein(nf) + "*" + ".nc")
    setfileoption("nc","SuppressClose",False) ; close too many files
    f      = addfiles(files(0:(nyear*12-1)),"r")
    ListSetType(f,"cat")
    ps   = f[:]->PS(:,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    dpi  = dpres_hybrid_ccm(ps,p0,hyai,hybi)  ;get dpi (time,lev,lat,lon)
    time = cd_calendar(f[:]->time,0)

    do nv = 0, nvar-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
    var_temp := f[:]->$var_name(nv)$(:,:,{lats:latn},{lonl:lonr})
    var = dim_sum_n_Wrap(var_temp*dpi,1) ;vertically integrated 

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
system("rm -f " + fileout(nf)+var_name(nv)+".nc")
ncdf = addfile(fileout(nf)+var_name(nv)+".nc","c")
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
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var_temp)
filevarattdef(ncdf,"season_ave",var_temp)
filevarattdef(ncdf,"month_vari",var_temp)
filevarattdef(ncdf,"season_vari",var_temp)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end do
end


