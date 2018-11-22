;******************************************************************************
;filename:  calc_Pb85_wave_flux_dzonal.ncl
;***********************************************************************
;According to Plumb1985-ncep-daily-10d-lpf.ncl 
;use seasonal climatology T,Z3 two avariables
;take zonal mean of every case as basic state
;calc every case monthly and seasonly climatology wave activity flux
;                                            by Ql_Ren
;                                           2018/11/08
;******************************************************************************
begin
;the path of the file
path     = "/users/yangsong3/renql/project/TP_NUDG"
case = (/"CTRL","NUDG24h","NUDG6h"/)
basic = 1
filein  = path + "/analysis/mdata/" + case + "-Clim_month_ave_" 
fileout = path + "/energy-20180417/mdata/" + case(basic) + "_Pl85_wave_flux_dzonal_seasonly.nc"
var_name = (/"season_ave"/) ;this is the var_name in data file
months   = (/"DJF","MAM","JJA","SON"/)
;months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)

;area: east Asia
lats = -90 ;-20
latn = 90  ;60
lonl = 0  ;20
lonr = 360 ;220
lev  = (/600,500,400,300,200/)

zf = addfile(filein(basic) + "Z3.nc","r")
tf = addfile(filein(basic) + "T.nc","r")
zvar = zf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})
tvar = tf->$var_name$(:,{lev},{lats:latn},{lonl:lonr})

ntime  = dimsizes(months)
nlat   = getfilevardimsizes(zf,"lat")
nlon   = getfilevardimsizes(zf,"lon")
nlev   = dimsizes(lev) 
lat   = zf->lat
lon   = zf->lon

;  Gas constant
gc=290
;  Gravitational acceleration
ga=9.80665

;  Radius of the earth
re=6378388

; scale height
sclhgt=8000.

; pi
pi = atan(1.0)*4.

; Coriolis parameter
f =  2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
f!0 = "lat"
f&lat = lat
f@_FillValue = -1.e+21

; missing for 10S - 10N
do ilat = 0, nlat-1
if (abs(lat(ilat) ).lt. 10. ) then
f(ilat)= f@_FillValue
end if
end do

; cosine
coslat = cos(lat(:)*pi/180.)

; zonal mean
tznl = dim_avg_n_Wrap(tvar(:,:,:,:),3)

; 1-D -> 4-D
leveltmp = conform_dims(dimsizes(tvar),lev,1)
coslattmp = conform_dims(dimsizes(tvar),coslat,2)
ftmp = conform_dims(dimsizes(tvar),f,2)

; 3-D -> 4-D
tznltmp = conform_dims( (/ntime,nlev,nlat,nlon/),tznl,(/0,1,2/))

; potential temperature
theta = tznltmp*(1000./leveltmp)^0.286

;K/m
dthetadz = center_finite_diff_n(theta,-sclhgt*log(lev/1000),False,0,1)

; Brunt Vaisala frequency
NN = (gc*(leveltmp/1000.)^0.286)/sclhgt * dthetadz

; deviations from zonal mean
psidev = dim_rmvmean(zvar)*ga/ftmp

;dpsidev/dlon
dpsidevdlon =  center_finite_diff_n(psidev,lon*pi/180.,True,0,3)

;ddpsidev/dlonlon
ddpsidevdlonlon =  center_finite_diff_n(dpsidevdlon,lon*pi/180.,True,0,3)

;dpsidev/dlat
dpsidevdlat = center_finite_diff_n(psidev, lat*pi/180., False,0,2)

;ddpsidev/dlonlat
ddpsidevdlonlat =  center_finite_diff_n(dpsidevdlon,lat*pi/180.,False,0,2)

;dpsidev/dz
dpsidevdz = center_finite_diff_n(psidev, -sclhgt*log(lev/1000),False,0,1)

;ddpsidev/dlondz
ddpsidevdlonz = center_finite_diff_n(dpsidevdlon, -sclhgt*log(lev/1000),False,0,1)

;x-component of (5.7)
Fx = leveltmp/1000./(2.*re*re*coslattmp)*(dpsidevdlon*dpsidevdlon-psidev*ddpsidevdlonlon)

;y-component of (5.7)
Fy = leveltmp/1000./(2.*re*re)*(dpsidevdlon*dpsidevdlat-psidev*ddpsidevdlonlat)

;z-component of (5.7)
Fz = leveltmp/1000.*ftmp*ftmp/(2.*NN*re)*(dpsidevdlon*dpsidevdz-psidev*ddpsidevdlonz)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "3-D wave-activity flux derived by Plumb (1985, JAS)" + \
                ", the basic state is zonal mean field of seasonal climatology" 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","lev","lat","lon"/)
dimSizes = (/ntime ,nlev ,nlat ,nlon/)
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

ncdf->month  = ispan(1,ntime,1)
ncdf->lev = (/lev/)
ncdf->lat = (/zvar&lat/)
ncdf->lon = (/zvar&lon/)
ncdf->Fx = (/Fx/)
ncdf->Fy = (/Fy/)
ncdf->Fz = (/Fz/)
ncdf->NN = (/NN/)
ncdf->psidev = (/psidev/)
end
