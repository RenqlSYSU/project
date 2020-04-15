;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","omega","hgt"/)

year_option = True;False ;
geo = False ;True ;if true, only can calc 10N 

fileout  = path2 + "ERA-Interim_DJF_zonal_momentum-year.nc"

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats =  10 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)

iyear = 0
itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
opt    = 0     ;used by center_finite_diff_n, no meanging 
loncyc = True  ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
nlev  = dimsizes(lev)
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
w     := new((/nyear,ntime,nlev,nlat,nlon/),float)
z     := new((/nyear,ntime,nlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end if
end do

;--------------------------------------
;DJF
;---------------------------------------------------
u2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
w2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
z2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    z2(ny,0:30,:,:,:) = z(ny,59:89,:,:,:)
    w2(ny,0:30,:,:,:) = w(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    z2(ny,31:89,:,:,:) = z(ny+1,0:58,:,:,:)
    w2(ny,31:89,:,:,:) = w(ny+1,0:58,:,:,:)
end do
u := u2
v := v2
z := z2
w := w2
printVarSummary(u)

;--------------------------------------
;Calculate seasonal averages and anomalies based on seasonal averages
;---------------------------------------------------
ca = 10.0 ;day
fca = 1.0/ca
fcf = fspan(0, 0.5, ntime/2) 
hgh_cf = ind(fcf.ge.fca)
low_cf = ind(fcf.lt.fca)

cf = ezfftf_n(u,itime)
term = cf 
term(:,:,hgh_cf,:,:,:) = 0.0  
low_u = ezfftb_n(term,0,itime+1)
term = cf 
term(:,:,low_cf,:,:,:) = 0.0  
hgh_u = ezfftb_n(term,0,itime+1)
printVarSummary(hgh_u)

cf = ezfftf_n(v,itime)
term = cf 
term(:,:,hgh_cf,:,:,:) = 0.0  
low_v = ezfftb_n(term,0,itime+1)
term = cf 
term(:,:,low_cf,:,:,:) = 0.0  
hgh_v = ezfftb_n(term,0,itime+1)
printVarSummary(hgh_v)

cf = ezfftf_n(w,itime)
term = cf 
term(:,:,hgh_cf,:,:,:) = 0.0  
low_w = ezfftb_n(term,0,itime+1)
term = cf 
term(:,:,low_cf,:,:,:) = 0.0  
hgh_w = ezfftb_n(term,0,itime+1)
printVarSummary(hgh_w)
delete([/cf,term,fcf,fca,hgh_cf,low_cf/])

coslat = conform(hgh_u,cos(lat),ilat)
tanlat = conform(hgh_u,tan(lat),ilat)

Al = dim_avg_n_Wrap(-center_finite_diff_n(low_u*low_u              ,lon,loncyc,opt,ilon)/a/coslat       ,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(low_v*low_u*coslat*coslat,lat,cyclic,opt,ilat)/a/coslat*coslat,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(low_w*low_u,lev*100,cyclic,opt,ilev),itime); + dim_avg_n_Wrap(low_u*low_v*tanlat/a,itime)
;Al  = dim_avg_n_Wrap(-(low_u/a/coslat)*center_finite_diff_n(low_u,lon,loncyc,opt,ilon),itime) + \
;      dim_avg_n_Wrap(-(low_v/a       )*center_finite_diff_n(low_u,lat,cyclic,opt,ilat),itime) + \
;      dim_avg_n_Wrap(-(low_w         )*center_finite_diff_n(low_u,lev*100,cyclic,opt,ilev),itime)
delete([/low_u,low_v,low_w/])
Ah = dim_avg_n_Wrap(-center_finite_diff_n(hgh_u*hgh_u              ,lon,loncyc,opt,ilon)/a/coslat       ,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(hgh_v*hgh_u*coslat*coslat,lat,cyclic,opt,ilat)/a/coslat*coslat,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(hgh_w*hgh_u,lev*100,cyclic,opt,ilev),itime); + dim_avg_n_Wrap(hgh_u*hgh_v*tanlat/a,itime)
;Ah  = dim_avg_n_Wrap(-(hgh_u/a/coslat)*center_finite_diff_n(hgh_u,lon,loncyc,opt,ilon),itime) + \
;      dim_avg_n_Wrap(-(hgh_v/a       )*center_finite_diff_n(hgh_u,lat,cyclic,opt,ilat),itime) + \
;      dim_avg_n_Wrap(-(hgh_w         )*center_finite_diff_n(hgh_u,lev*100,cyclic,opt,ilev),itime)
delete([/hgh_u,hgh_v,hgh_w/])

;--------------------------------------
;Calculate seasonal averages and anomalies based on seasonal averages
;---------------------------------------------------
clm_u = dim_avg_n(u,itime)
clm_v = dim_avg_n(v,itime) 
clm_z = dim_avg_n(z,itime) 
clm_w = dim_avg_n(w,itime) 
edd_u = u - conform(u, clm_u, (/0,2,3,4/))
edd_v = v - conform(v, clm_v, (/0,2,3,4/))   
edd_w = w - conform(w, clm_w, (/0,2,3,4/))  
printVarSummary(clm_u)
delete([/u,v,z,w/])

;==========================================================
;calc forcing factors
;============================================================
;A  = dim_avg_n_Wrap(-(edd_u/a/coslat)*center_finite_diff_n(edd_u,lon,loncyc,opt,ilon),itime) + \
;     dim_avg_n_Wrap(-(edd_v/a       )*center_finite_diff_n(edd_u,lat,cyclic,opt,ilat),itime) + \
;     dim_avg_n_Wrap(-(edd_w         )*center_finite_diff_n(edd_u,lev*100,cyclic,opt,ilev),itime)
A  = dim_avg_n_Wrap(-center_finite_diff_n(edd_u*edd_u              ,lon,loncyc,opt,ilon)/a/coslat       ,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(edd_v*edd_u*coslat*coslat,lat,cyclic,opt,ilat)/a/coslat*coslat,itime) + \
     dim_avg_n_Wrap(-center_finite_diff_n(edd_w*edd_u,lev*100,cyclic,opt,ilev),itime); + dim_avg_n_Wrap(edd_u*edd_v*tanlat/a,itime)
printVarSummary(A)
print(A(0,5,30,30))
delete([/edd_u,edd_v,edd_w/])
delete([/coslat,tanlat/])

coslat4d = conform(clm_u,cos(lat),ilat-1)
tanlat4d = conform(clm_u,tan(lat),ilat-1)
f0  = conform(clm_u,2*(2*pi/24.0/3600.0)*sin(lat),ilat-1)
fva = f0*clm_v - center_finite_diff_n(clm_z,lon,loncyc,opt,ilon-1)/a/coslat4d
delete([/clm_z,f0/])

adv  = -(clm_u/a/coslat4d)*center_finite_diff_n(clm_u,lon,loncyc,opt,ilon-1)\
       -(clm_v/a         )*center_finite_diff_n(clm_u,lat,cyclic,opt,ilat-1)\
       -(clm_w           )*center_finite_diff_n(clm_u,lev*100,cyclic,opt,ilev-1)
printVarSummary(adv)
print(adv(0,5,30,30))

cur = clm_u*clm_v*tanlat4d/a - clm_u*clm_w/a
delete([/clm_u,clm_v,clm_w/])
delete([/coslat4d,tanlat4d/])

Ah@long_name = "convergence of zonal-momentum fuxes by synoptic eddies" 
Ah@units     = "m*s^-2" 
Al@long_name = "convergence of zonal-momentum fuxes by low-frequency eddies" 
Al@units     = "m*s^-2" 
A@long_name = "convergence of zonal-momentum fuxes by TE" 
A@units     = "m*s^-2" 
fva@long_name = "ageostrophic acceleration associated with meridional overturning circulations" 
fva@units     = "m*s^-2" 
adv@long_name = "advection by the time-mean fow" 
adv@units     = "m*s^-2" 
cur@long_name = "curvature" 
cur@units     = "m*s^-2" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "DJF mean three forcing factors, Qd, Qeddy, A, fixed f0, geo"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year" ,"lev","lat","lon"/);
dimSizes = (/nyear-1,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"Ah"  ,typeof(Ah),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Al"  ,typeof(Al),(/"year","lev","lat","lon"/))
filevardef(ncdf,"A"   ,typeof(A),(/"year","lev","lat","lon"/))
filevardef(ncdf,"fva" ,typeof(fva),(/"year","lev","lat","lon"/))
filevardef(ncdf,"adv" ,typeof(adv),(/"year","lev","lat","lon"/))
filevardef(ncdf,"cur" ,typeof(cur),(/"year","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"Ah" ,Ah  )
filevarattdef(ncdf,"Al" ,Al  )
filevarattdef(ncdf,"A"  ,A   )
filevarattdef(ncdf,"fva",fva )
filevarattdef(ncdf,"adv",adv )
filevarattdef(ncdf,"cur",cur )

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->Ah   = (/Ah/)
ncdf->Al   = (/Al/)
ncdf->A    = (/A /)
ncdf->fva  = (/fva/)
ncdf->adv  = (/adv/)
ncdf->cur  = (/cur/)
end

