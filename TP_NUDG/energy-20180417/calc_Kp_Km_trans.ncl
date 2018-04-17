;***************************************************************************
;filename:      calc_Kp_Km_trans.ncl
;use the initial data to calculate the transfer between perturbation kinetic energy
;and mean flow kinetic energy which is referenced in Dynamic Meteorology P210.
;The ncl functions have been used are: dpres_hybrid_ccm,dim_sum_n_Wrap,dim_avg_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = (/"/users/yangsong3/L_Zealot/F/AMIP_C5PM/exe/AMIP_C5PM.cam.h1.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/pre_data_24h/AMIP_C5PM_TP_NUDG.cam.h1.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/AMIP_C5PM_TP_NUDG.cam.h1."/)
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Kp_Km_trans_10avg_vert_int"
case    = (/"CTRL","NUDG24h","NUDG6h"/)
years   = ispan(1979,2005,1)   ;the calculate time
ncase = dimsizes(case)

;area of field variables,larger then the plan
lats = -20
latn = 60
lonl = 23 ;30
lonr = 157 ;150

files1 = systemfunc("ls " + filein(0) + years(0) + "*" + ".nc")
f1     = addfile(files1(0),"r")
hyai   = f1->hyai  ;hyai(ilev),ilev = lev+1,hybrid A coefficient at layer interfaces
hybi   = f1->hybi
hyam   = f1->hyam  ;hybm(lev), hybrid A coefficient at layer midpoints
hybm   = f1->hybm
p0     = 100000   ;a scalar ,double,reference pressure,unit is Pa
printVarSummary(hyai)
printVarSummary(hyam)

do nf = 0,ncase-1,1 
;==========================================================
;read the data
;============================================================
files  = systemfunc("ls " + filein(nf) + years(0) + "*" + ".nc")
f      = addfiles(files(0:2),"r")
ListSetType(f,"cat")
uwnd   = f[:]->U(:,:,{lats:latn},{lonl:lonr})
vwnd   = f[:]->V(:,:,{lats:latn},{lonl:lonr})
wwnd   = f[:]->OMEGA(:,:,{lats:latn},{lonl:lonr})
ps     = f[:]->PS(:,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
printVarSummary(ps)
printVarSummary(uwnd)

nlat  = dimsizes(uwnd&lat)
nlon  = dimsizes(uwnd&lon)
nlev  = dimsizes(uwnd&lev)
ntime = dimsizes(uwnd&time)
var1  = new((/ntime,nlev,nlat,nlon/),double)
var2  = new((/ntime,nlev,nlat,nlon/),double)
var3  = new((/ntime,nlev,nlat,nlon/),double)
var4  = new((/ntime,nlev,nlat,nlon/),double)

;calculates the pressure differences of a hybrid coordinate system.
dpi  = dpres_hybrid_ccm(ps,p0,hyai,hybi)  ;get dpi (time,lev,lat,lon)
plev = pres_hybrid_ccm( ps,p0,hyam,hybm)  ;get plev(time,lev,lat,lon),used to calc the central difference
printVarSummary(dpi)
printVarSummary(plev)

;calculat zonal mean
zm_u = dim_avg_n_Wrap(uwnd,3) ;get zm_u(time,lev,lat) 
zm_v = dim_avg_n_Wrap(vwnd,3)  
zm_w = dim_avg_n_Wrap(wwnd,3)  
zm_plev = dim_avg_n_Wrap(plev,3)  

;calculate perturbation
zp_u = dim_rmvmean_n_Wrap(uwnd,3) ;get zp_u(time,lev,lat,lon)
zp_v = dim_rmvmean_n_Wrap(vwnd,3) 
zp_w = dim_rmvmean_n_Wrap(wwnd,3) 

var1  = 0
var2  = 0
var3  = 0
var4  = 0
do nlo = 5, nlon-6, 1
    var1(:,:,:,nlo) = dim_avg_n_Wrap(zp_u(:,:,:,nlo-5:nlo+5)*zp_v(:,:,:,nlo-5:nlo+5),3) \
                      *center_finite_diff_n(zm_u, uwnd&lat, False, 0, 2) ;center_finite_diff_n (data, coordinates, rcyclic, opt, dim)
    var2(:,:,:,nlo) = dim_avg_n_Wrap(zp_u(:,:,:,nlo-5:nlo+5)*zp_w(:,:,:,nlo-5:nlo+5),3) \
                       *center_finite_diff_n(zm_u, zm_plev, False, 0, 1) ;center_finite_diff_n (data, coordinates, rcyclic, opt, dim)
    var3(:,:,:,nlo) = dim_avg_n_Wrap(zp_v(:,:,:,nlo-5:nlo+5)*zp_v(:,:,:,nlo-5:nlo+5),3) \
                       *center_finite_diff_n(zm_v, uwnd&lat, False, 0, 2) ;center_finite_diff_n (data, coordinates, rcyclic, opt, dim)
    var4(:,:,:,nlo) = dim_avg_n_Wrap(zp_v(:,:,:,nlo-5:nlo+5)*zp_w(:,:,:,nlo-5:nlo+5),3) \
                       *center_finite_diff_n(zm_v, zm_plev, False, 0, 1) ;center_finite_diff_n (data, coordinates, rcyclic, opt, dim)
end do
kpkm =dim_sum_n_Wrap(( var1+var2+var3+var4)*dpi,1) ;Vertical integration,get(ntime,nlat,nlon)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+years(0) + ".nc")
ncdf = addfile(fileout(nf)+years(0) + ".nc","c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "the transfer between perturbation kinetic energyc and mean flow kinetic energy "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lat" ,"lon"/)
dimSizes = (/ntime ,nlat  ,nlon /)
dimUnlim = (/True  ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(uwnd&time),(/"time"/))
filevardef(ncdf,"lat" ,typeof(uwnd&lat),getvardims(uwnd&lat))
filevardef(ncdf,"lon" ,typeof(uwnd&lon),getvardims(uwnd&lon))
filevardef(ncdf,"kpkm",typeof(kpkm),(/"time","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",uwnd&lat)
filevarattdef(ncdf,"lat" ,uwnd&lat)
filevarattdef(ncdf,"lon" ,uwnd&lon)
kpkm@long_name = "vertical integrated and ten grids average, >0 mean Kp convert to Km"
kpkm@units     = "J"
filevarattdef(ncdf,"kpkm",kpkm)

ncdf->time = (/uwnd&time/)
ncdf->lat  = (/uwnd&lat/)
ncdf->lon  = (/uwnd&lon/)
ncdf->kpkm = (/kpkm/)
end do
end

