;******************************************************************************
;filename:  calc_TN01_wave_flux_dclim.ncl
;According to TakayaNakamura2001-ncep-monthly.ncl file
;use Z3 and climatology U,V,T,Z3, take climatology of NUDG case as basic state
;calc CTRL monthly and seasonly climatology wave activity flux
;                                            by Ql_Ren
;                                           2018/11/07
;modified:
;2018/11/08: calc every year daily anomaly Z3 and wave activity flux
;            take each case's climatology as basic state
;            then calc seasonal climatology wave activity flux
;******************************************************************************
begin
;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
case = (/"CTRL","NUDG24h","NUDG6h"/)
basic = 1
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout  = path + "/energy-20180417/mdata/" + case(basic) + "_TN01_wave_flux_dclim_seasonly.nc"
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = 10 ;-90 ;-20
latn = 90  ;60
lonl = 0  ;20
lonr = 360 ;220
lev  = (/600,500,400,300,200/)

zf = addfile(filein(basic) + years(0) + ".daily.Z3.nc","r")
var = zf->Z3(0,{lev(0)},{lats:latn},{lonl:lonr})
nyear  = dimsizes(years)
ntime  = 365
nlat   = dimsizes(var&lat) 
nlon   = dimsizes(var&lon)
nlev   = dimsizes(lev) 
lat   = var&lat
lon   = var&lon
delete(var)
zvar   = new((/nyear,365,nlev,nlat,nlon/),float)
tvar   = new((/nyear,365,nlev,nlat,nlon/),float)
uvar   = new((/nyear,365,nlev,nlat,nlon/),float)
vvar   = new((/nyear,365,nlev,nlat,nlon/),float)
dlev = 1
dlat = 2
dlon = 3

do nt = 0, nyear-1, 1
zf = addfile(filein(basic) + years(nt) + ".daily.Z3.nc","r")
tf = addfile(filein(basic) + years(nt) + ".daily.T.nc","r")
uf = addfile(filein(basic) + years(nt) + ".daily.U.nc","r")
vf = addfile(filein(basic) + years(nt) + ".daily.V.nc","r")
zvar(nt,:,:,:,:) = zf->Z3(:,{lev},{lats:latn},{lonl:lonr})
tvar(nt,:,:,:,:) = tf->T(:,{lev},{lats:latn},{lonl:lonr})
uvar(nt,:,:,:,:) = uf->U(:,{lev},{lats:latn},{lonl:lonr})
vvar(nt,:,:,:,:) = vf->V(:,{lev},{lats:latn},{lonl:lonr})
end do
delete(tf)
delete(uf)
delete(vf)
printVarSummary(zvar)

bzvar = dim_avg_n_Wrap(zvar,0)
btvar = dim_avg_n_Wrap(tvar,0)
buvar = dim_avg_n_Wrap(uvar,0)
bvvar = dim_avg_n_Wrap(vvar,0)
delete(uvar)
delete(tvar)
delete(vvar)
printVarSummary(buvar)

Fx     = new((/nyear,365,nlev,nlat,nlon/),double);float)
Fy     = new((/nyear,365,nlev,nlat,nlon/),double);float)
Fz     = new((/nyear,365,nlev,nlat,nlon/),double);float)
psidev = new((/nyear,365,nlev,nlat,nlon/),double);float)

;  Gas constant,unit is J/Kg/K
gc=290 
;  Gravitational acceleration
ga=9.80665
;  Radius of the earth, unit is m
re=6378388
; scale height
sclhgt=8000.
; pi
pi = atan(1.0)*4.

; Coriolis parameter
f = tofloat( 2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:)))
f!0 = "lat"
f&lat = lat
f@_FillValue = zvar@_FillValue

; missing for 10S - 10N
;do ilat = 0, nlat-1
;if (abs(lat(ilat) ).lt. 10. ) then
;f(ilat)= f@_FillValue
;end if
;end do

; cosine
coslat = cos(lat(:)*pi/180.)

; 1-D -> 4-D
leveltmp = conform_dims(dimsizes(bzvar),lev,dlev)
coslattmp = conform_dims(dimsizes(bzvar),coslat,dlat)
ftmp = conform_dims(dimsizes(bzvar),f,dlat)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(btvar*(1000./leveltmp)^0.286,-sclhgt*log(lev/1000),False,0,dlev)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = btvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(buvar^2 + bvvar^2)
cumag@_FillValue = buvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

do nt = 0, nyear-1, 1
;height anomaly
zavar = zvar(nt,:,:,:,:) - bzvar
printVarSummary(zavar)
printVarSummary(ga)
printVarSummary(ftmp)

; QG steam function for anomaly
psidev(nt,:,:,:,:) = zavar*ga /ftmp

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev(nt,:,:,:,:),lon*pi/180.,True,0,dlon)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,dlon)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev(nt,:,:,:,:), lat*pi/180., False,0,dlat)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,dlat)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,dlat)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev(nt,:,:,:,:), -sclhgt*log(lev/1000),False,0,dlev)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(lev/1000),False,0,dlev)

;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(lev/1000),False,0,dlev)

xuterm = (dpsidevdlon*dpsidevdlon - psidev(nt,:,:,:,:)*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev(nt,:,:,:,:)*ddpsidevdlonlat)
;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev(nt,:,:,:,:)*ddpsidevdlatlat)

zuterm = (dpsidevdlon*dpsidevdz - psidev(nt,:,:,:,:)*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev(nt,:,:,:,:)*ddpsidevdlatz)

; Mask out where westerlies is small or negative (less than 5 m/s).
Fx(nt,:,:,:,:) = mask( leveltmp/1000./(2.*cumag*re*re)*( buvar/coslattmp * xuterm + bvvar * xvterm), buvar.lt.5,False)
Fy(nt,:,:,:,:) = mask(leveltmp/1000./(2.*cumag*re*re)*( buvar*xvterm + coslattmp*bvvar*yvterm),  buvar.lt.5,False)
Fz(nt,:,:,:,:) = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(buvar*zuterm + coslattmp*bvvar*zvterm), buvar.lt.5, False)
end do
printVarSummary(Fx)
copy_VarMeta(zvar,Fx    )
copy_VarMeta(zvar,Fy    )
copy_VarMeta(zvar,Fz    )
copy_VarMeta(zvar,psidev)
delete(leveltmp) 
delete(coslattmp)
delete(ftmp)
delete(cumag)
delete(zavar)

;==============================================================================
;calc wave activity flux climatology
;==============================================================================
month_Fx     = new((/nyear,12,nlev,nlat,nlon/),double);float)
month_Fy     = new((/nyear,12,nlev,nlat,nlon/),double);float)
month_Fz     = new((/nyear,12,nlev,nlat,nlon/),double);float)
month_psidev = new((/nyear,12,nlev,nlat,nlon/),double);float)
season_Fx     = new((/nyear,4,nlev,nlat,nlon/),double);float)
season_Fy     = new((/nyear,4,nlev,nlat,nlon/),double);float)
season_Fz     = new((/nyear,4,nlev,nlat,nlon/),double);float)
season_psidev = new((/nyear,4,nlev,nlat,nlon/),double);float)

time = cd_calendar(zf->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    month_Fx(:,nm,:,:,:) = dim_avg_n_Wrap(Fx(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    month_Fy(:,nm,:,:,:) = dim_avg_n_Wrap(Fy(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    month_Fz(:,nm,:,:,:) = dim_avg_n_Wrap(Fz(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    month_psidev(:,nm,:,:,:) = dim_avg_n_Wrap(psidev(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

;----------calculate the season average-----------------------
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    select_time := ind((time(:,1).eq.nm+1).or.(time(:,1).eq.nm).or.(time(:,1).eq.11))
    else
    select_time := ind((time(:,1).eq.nm+1).or.(time(:,1).eq.nm).or.(time(:,1).eq.nm-1))
    end if
    season_Fx(:,ns,:,:,:)     = dim_avg_n_Wrap(Fx(:,select_time,:,:,:),1) ;get (nyear,nseason,nlat,nlon)
    season_Fy(:,ns,:,:,:)     = dim_avg_n_Wrap(Fy(:,select_time,:,:,:),1) ;get (nyear,nseason,nlat,nlon)
    season_Fz(:,ns,:,:,:)     = dim_avg_n_Wrap(Fz(:,select_time,:,:,:),1) ;get (nyear,nseason,nlat,nlon)
    season_psidev(:,ns,:,:,:) = dim_avg_n_Wrap(psidev(:,select_time,:,:,:),1) ;get (nyear,nseason,nlat,nlon)
end do
delete(Fx)
delete(Fy)

month_ave_Fx      = dim_avg_n_Wrap(month_Fx    ,0) ;get (nmonth,nlat,nlon)
month_ave_Fy      = dim_avg_n_Wrap(month_Fy    ,0) ;get (nmonth,nlat,nlon)
month_ave_Fz      = dim_avg_n_Wrap(month_Fz    ,0) ;get (nmonth,nlat,nlon)
month_ave_psidev  = dim_avg_n_Wrap(month_psidev,0) ;get (nmonth,nlat,nlon)

month_vari_Fx     = dim_variance_n_Wrap(month_Fx    ,0) ;get (nmonth,nlat,nlon)
month_vari_Fy     = dim_variance_n_Wrap(month_Fy    ,0) ;get (nmonth,nlat,nlon)
month_vari_Fz     = dim_variance_n_Wrap(month_Fz    ,0) ;get (nmonth,nlat,nlon)
month_vari_psidev = dim_variance_n_Wrap(month_psidev,0) ;get (nmonth,nlat,nlon)

season_ave_Fx      = dim_avg_n_Wrap(season_Fx    ,0) ;get (nseason,nlat,nlon)
season_ave_Fy      = dim_avg_n_Wrap(season_Fy    ,0) ;get (nseason,nlat,nlon)
season_ave_Fz      = dim_avg_n_Wrap(season_Fz    ,0) ;get (nseason,nlat,nlon)
season_ave_psidev  = dim_avg_n_Wrap(season_psidev,0) ;get (nseason,nlat,nlon)

season_vari_Fx     = dim_variance_n_Wrap(season_Fx    ,0) ;get (nseason,nlat,nlon)
season_vari_Fy     = dim_variance_n_Wrap(season_Fy    ,0) ;get (nseason,nlat,nlon)
season_vari_Fz     = dim_variance_n_Wrap(season_Fz    ,0) ;get (nseason,nlat,nlon)
season_vari_psidev = dim_variance_n_Wrap(season_psidev,0) ;get (nseason,nlat,nlon)

printVarSummary(month_ave_Fx)
printVarSummary(month_vari_Fx)
printVarSummary(month_ave_Fz)
printVarSummary(month_vari_Fz)

delete(month_Fx     )
delete(month_Fy     )
delete(month_Fz     )
delete(month_psidev )
delete(season_Fx    )
delete(season_Fy    )
delete(season_Fz    )
delete(season_psidev)
;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)" + \
               ", the basic state is itself daily climatology"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","season","lev","lat","lon"/)
dimSizes = (/  12   ,  4     ,nlev ,nlat ,nlon/)
dimUnlim = (/False  ,False   ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(zvar&lat),getvardims(zvar&lat))
filevardef(ncdf,"lon",typeof(zvar&lon),getvardims(zvar&lon))
filevardef(ncdf,"month_ave_Fx"     ,typeof(month_ave_Fx     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_ave_Fy"     ,typeof(month_ave_Fy     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_ave_Fz"     ,typeof(month_ave_Fz     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_ave_psidev" ,typeof(month_ave_psidev ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave_Fx"    ,typeof(season_ave_Fx    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_ave_Fy"    ,typeof(season_ave_Fy    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_ave_Fz"    ,typeof(season_ave_Fz    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_ave_psidev",typeof(season_ave_psidev),(/"season","lev","lat","lon"/))
filevardef(ncdf,"month_vari_Fx"     ,typeof(month_vari_Fx     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari_Fy"     ,typeof(month_vari_Fy     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari_Fz"     ,typeof(month_vari_Fz     ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari_psidev" ,typeof(month_vari_psidev ) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_vari_Fx"    ,typeof(season_vari_Fx    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_vari_Fy"    ,typeof(season_vari_Fy    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_vari_Fz"    ,typeof(season_vari_Fz    ),(/"season","lev","lat","lon"/))
filevardef(ncdf,"season_vari_psidev",typeof(season_vari_psidev),(/"season","lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",zvar&level)
filevarattdef(ncdf,"lat",zvar&lat)
filevarattdef(ncdf,"lon",zvar&lon)
Fz@units = "m^2/s^2"
psidev@units = "m^2/s"
filevarattdef(ncdf,"month_ave_Fx"     , Fz )
filevarattdef(ncdf,"month_ave_Fy"     , Fz )
filevarattdef(ncdf,"month_ave_Fz"     , Fz )
filevarattdef(ncdf,"month_ave_psidev" ,psidev)
filevarattdef(ncdf,"season_ave_Fx"    , Fz )
filevarattdef(ncdf,"season_ave_Fy"    , Fz )
filevarattdef(ncdf,"season_ave_Fz"    , Fz )
filevarattdef(ncdf,"season_ave_psidev",psidev)
filevarattdef(ncdf,"month_vari_Fx"     , Fz )
filevarattdef(ncdf,"month_vari_Fy"     , Fz )
filevarattdef(ncdf,"month_vari_Fz"     , Fz )
filevarattdef(ncdf,"month_vari_psidev" ,psidev)
filevarattdef(ncdf,"season_vari_Fx"    , Fz )
filevarattdef(ncdf,"season_vari_Fy"    , Fz )
filevarattdef(ncdf,"season_vari_Fz"    , Fz )
filevarattdef(ncdf,"season_vari_psidev",psidev)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lev = (/lev/)
ncdf->lat = (/zvar&lat/)
ncdf->lon = (/zvar&lon/)
ncdf->month_ave_Fx      = (/month_ave_Fx     /)
ncdf->month_ave_Fy      = (/month_ave_Fy     /)
ncdf->month_ave_Fz      = (/month_ave_Fz     /)
ncdf->month_ave_psidev  = (/month_ave_psidev /)
ncdf->season_ave_Fx     = (/season_ave_Fx    /)
ncdf->season_ave_Fy     = (/season_ave_Fy    /)
ncdf->season_ave_Fz     = (/season_ave_Fz    /)
ncdf->season_ave_psidev = (/season_ave_psidev/)
ncdf->month_vari_Fx      = (/month_vari_Fx     /)
ncdf->month_vari_Fy      = (/month_vari_Fy     /)
ncdf->month_vari_Fz      = (/month_vari_Fz     /)
ncdf->month_vari_psidev  = (/month_vari_psidev /)
ncdf->season_vari_Fx     = (/season_vari_Fx    /)
ncdf->season_vari_Fy     = (/season_vari_Fy    /)
ncdf->season_vari_Fz     = (/season_vari_Fz    /)
ncdf->season_vari_psidev = (/season_vari_psidev/)
end
