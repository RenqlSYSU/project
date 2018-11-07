;******************************************************************************
;filename:  calc_TN01_wave_flux_dclim.ncl
;According to TakayaNakamura2001-ncep-monthly.ncl file
;use U,V,T,Z3 four avariables, take climatology of NUDG case as basic state
;calc CTRL monthly and seasonly climatology wave activity flux
;                                            by Ql_Ren
;                                           2018/11/08
;******************************************************************************
begin
;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
;path     = "/home/ys17-19/renql/project/TP_NUDG"
case = (/"CTRL","NUDG24h","NUDG6h"/)
basic = 2
filein  = path + "/analysis/mdata/" + case + "-Clim_month_ave_" 
fileout = path + "/energy-20180417/mdata/CTRL_b" + case(basic) + "_TN01_wave_flux_dclim_seasonly.nc"
var_name = (/"season_ave"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;area: east Asia
lats = -90 ;-20
latn = 90  ;60
lonl = 0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600,500,400,300,200/)

zf = addfile(filein(0) + "Z3.nc","r")
zvar = zf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})

;read basic state
bzf = addfile(filein(basic) + "Z3.nc","r")
btf = addfile(filein(basic) + "T.nc","r")
buf = addfile(filein(basic) + "U.nc","r")
bvf = addfile(filein(basic) + "V.nc","r")
bzvar = bzf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})
btvar = bzf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})
buvar = bzf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})
bvvar = bzf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})

nmonth = dimsizes(months)
nlat   = getfilevardimsizes(zf,"lat")
nlon   = getfilevardimsizes(zf,"lon")
nlev    = getfilevardimsizes(zf,"lev")
lat   = zf->lat
lon   = zf->lon

;height anomaly
zavar = zvar - bzvar

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
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
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
leveltmp = conform_dims(dimsizes(zvar),lev,1)
coslattmp = conform_dims(dimsizes(zvar),coslat,2)
ftmp = conform_dims(dimsizes(zvar),f,2)

; vertical gradient of potential temperature  (K/m)
dthetadz = center_finite_diff_n(btvar*(1000./leveltmp)^0.286,-sclhgt*log(lev/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
NN@_FillValue = btvar@_FillValue
NN = where(NN .gt. 0, NN, NN@_FillValue)

; magnitude of climatological wind
cumag = sqrt(buvar^2 + bvvar^2)
cumag@_FillValue = buvar@_FillValue
cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

; QG steam function for anomaly
psidev = zavar*ga /ftmp

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;ddpsidev/dlatdlat
ddpsidevdlatlat = center_finite_diff_n(dpsidevdlat, lat*pi/180.,False,0,2)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(lev/1000),False,0,1)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(lev/1000),False,0,1)

;ddpsidev/dlatdz
ddpsidevdlatz = center_finite_diff_n(dpsidevdlat, -sclhgt*log(lev/1000),False,0,1)

xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
;yuterm = xvterm
yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

; Mask out where westerlies is small or negative (less than 5 m/s).
;x-component of (38)
Fx = mask( leveltmp/1000./(2.*cumag*re*re)*( buvar/coslattmp * xuterm + bvvar * xvterm), buvar.lt.5,False)
;y-component 
Fy = mask(leveltmp/1000./(2.*cumag*re*re)*( buvar*xvterm + coslattmp*bvvar*yvterm),  buvar.lt.5,False)
;z-component 
Fz = mask(leveltmp/1000.*ftmp*ftmp/(2.*cumag*NN*re)*(buvar*zuterm + coslattmp*bvvar*zvterm), buvar.lt.5, False)

printVarSummary(Fx)

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
                ", the basic state is NUDG seasonly climatology"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","lev","lat","lon"/)
dimSizes = (/nmonth ,nlev ,nlat ,nlon/)
dimUnlim = (/False  ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(zvar&lat),getvardims(zvar&lat))
filevardef(ncdf,"lon",typeof(zvar&lon),getvardims(zvar&lon))
filevardef(ncdf,"Fx" ,typeof(Fx) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"Fy" ,typeof(Fy) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"Fz" ,typeof(Fz) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"NN" ,typeof(NN) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"psidev" ,typeof(psidev) ,(/"month" ,"lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",zvar&lev)
filevarattdef(ncdf,"lat",zvar&lat)
filevarattdef(ncdf,"lon",zvar&lon)
Fx@units = "m^2/s^2"
Fx@units = "m^2/s^2"
Fz@units = "m^2/s^2"
psidev@units = "m^2/s"
NN@units = "1/s^2"
filevarattdef(ncdf,   "Fx"  , Fx )
filevarattdef(ncdf,   "Fy"  , Fy )
filevarattdef(ncdf,   "Fz"  , Fz )
filevarattdef(ncdf,   "NN"  , NN )
filevarattdef(ncdf,"psidev" ,psidev)

ncdf->month  = ispan(1,nmonth,1)
ncdf->lev = (/lev/)
ncdf->lat = (/zvar&lat/)
ncdf->lon = (/zvar&lon/)
ncdf->Fx = (/Fx/)
ncdf->Fy = (/Fy/)
ncdf->Fz = (/Fz/)
ncdf->NN = (/NN/)
ncdf->psidev = (/psidev/)
end
