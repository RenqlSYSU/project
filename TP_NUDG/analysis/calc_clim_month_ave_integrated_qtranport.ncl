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
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_ave_"
var_name  = (/"U","V","Q","PS"/)
calc_name = (/"qu","qv","qdivg"/)
years    = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/1000,925,850,700,500/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(vars&level)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var_u      = new((/nyear,365,nlev,nlat,nlon/),float)
var_v      = new((/nyear,365,nlev,nlat,nlon/),float)
var_q      = new((/nyear,365,nlev,nlat,nlon/),float)
var_ps     = new((/nyear,365,nlat,nlon/),float)
var        = new((/nyear,365,nlat,nlon/),float)
var_month  = new((/nyear,12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)

do nf = 0 , ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
nv = 0
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_u(nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years

nv = 1
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_v(nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years

nv = 2
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_q(nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years
var_q = var_q*1000  ;convert unit from kg/kg to g/kg

nv =3
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_ps(nt,:,:,:) = f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years
var_ps = var_ps/100 ;convert unit from Pa to hPa

var_u = var_u*var_q/g
var_v = var_v*var_q/g
var_divg = uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 2) ;calculate the divergence

do nvc = 2, dimsizes(calc_name)-1,1
if (nvc.eq.0) then
    var = wgt_vert_avg_beta ( var_u&level, var_u, var_ps, 0, 0 )
    ;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum
    var@long_name = "Total (vertically integrated) zonal water vapor transmission"
    var@units     = "g/(s*m)"
else if (nvc.eq.1) then
    var = wgt_vert_avg_beta ( var_v&level, var_v, var_ps, 0, 0 )
    ;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum
    var@long_name = "Total (vertically integrated) meridional water vapor transmission"
    var@units     = "g/(s*m)"
else
    var = wgt_vert_avg_beta ( var_u&level, var_divg, var_ps, 0, 0 )
    ;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum
    var@long_name = "Total (vertically integrated) divergence of moisture flux"
    var@units     = "g/(s*m*m)"
end if
end if

;----------calculate the month and season average-----------------------
time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

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
system("rm -f " + fileout(nf)+calc_name(nvc)+".nc")
ncdf = addfile(fileout(nf)+calc_name(nvc)+".nc","c")
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
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

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


