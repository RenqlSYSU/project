;***************************************************************************
;filename:      calc_season_day_vari_preci.ncl
;output_file:   Season_day_vari_preci.nc 
;read the 27 years daily U,V,T at East Asia and then
;get an array of (nvar,nyear,nday,nlevel,nlat,nlon),then calculate the variance in
;the dim of nday to get (nvar,nyear,nmonth,nlevel,nlat,nlon),then calculate
;the average in the dim of nyear.
;                                            by Ql_Ren
;                                           2018/01/17
;******************************************************************************
begin
;the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Season_day_vari_rmvmean_"
var_name = (/"TMQ","LHFLX","SHFLX"/)
case     = (/"AMIP_CTRL","TP_NUDG24h","TP_NUDG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nyear= dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nday  = 365
nmon  = 12
nsea  = 4
var       = new((/nyear,nday,nlat,nlon/),float)
vari_temp = new((/nyear,nsea,nlat,nlon/),float)
vari      = new((/      nsea,nlat,nlon/),float)

;do nv = 0 ,nvar-1,1 
do nv = 0 ,0,1 
do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nt,:,:,:) = f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
end do  ;the loop of the 27 years
var = dim_rmvmean_n_Wrap(var,0)

;-------calc the daily variance in 12 month-------------
time = cd_calendar(f->time,0)
do ns = 0,nsea-1,1
    if(ns.eq.0)
    select_time := ind(time(:,1).lt.3.or.time(:,1).eq.12)
    else
    nm_s = 3*ns
    nm_e = 3*(ns+1)-1
    select_time := ind(time(:,1).le.nm_e.and.time(:,1).ge.nm_s)
    end if
    ntime = dimsizes(select_time)
    print(ntime)
    vari_temp(:,ns,:,:) = dim_variance_n_Wrap(var(:,select_time,:,:),1) ;get (nvar,nyear,nmonth,nlev,nlat,nlon)
end do
vari = dim_avg_n_Wrap(vari_temp,0)  ;get (nvar,nmonth,nlev,nlat,nlon)
printVarSummary(vari)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc)+var_name(nv)+".nc")
ncdf = addfile(fileout(nc)+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "daily variance in each season, averaged over 27 years"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"season","lat","lon"/)
dimSizes = (/nsea    ,nlat ,nlon/)
dimUnlim = (/False   ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"season_vari",typeof(vari),(/"season","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
vari@long_name = "daily variance of " + var_name(nv) + " in 4 seasons"
vari@units   = ""
filevarattdef(ncdf,"season_vari",vari)

ncdf->season  = ispan(1,4,1)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->season_vari = (/vari/)
end do
end do
end

