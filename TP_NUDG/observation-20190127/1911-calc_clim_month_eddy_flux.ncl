;***************************************************************************
;use the post-processing data to calculate Qd, Qeddy and A
;1. calc climate and transient U V T W
;2. calc transient Qeddy and A, climate Qd 
;3. calc climatology daily Qeddy and A
;4. calc climatology daily forcing induced by Qeddy and Qd,
;   which can influence geopotential tendency
;
;The ncl functions have been used are: 
;   dim_avg_n_Wrap, center_finite_diff_n, conform
;   dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/01/25
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","omega","air"/)
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter

;parameter for bw_bandpass_filter
ca = 7.0 ;day
cb = 2.5 

if(option.eq.0) then
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_eddy_flux.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_eddy_flux-" + cb + "-" + ca + "day.nc"
end if

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats = -2 
latn = 90 
lonl = 0  
lonr = 360
lev  = (/1000,975,950,925,900, 875,850,825,800,775, 750,700,650,600,550, 500,450,400,350,300, 250,225,200,175,150, 125,100/)
month_s = 1
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
nmonth= 4
Qeddy_m  = new((/ 3,nmonth,nlev,nlat,nlon/),float)
A_m      = new((/ 2,nmonth,nlev,nlat,nlon/),float)
Tflux_m  = new((/ 3,nmonth,nlev,nlat,nlon/),float)
Sflux_m  = new((/ 3,nmonth,nlev,nlat,nlon/),float)
Mflux_m  = new((/ 3,nmonth,nlev,nlat,nlon/),float)

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
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
do nm = 0, nmonth-1, 1
month_s = 2*nm + 1
month_e = 2*nm + 2
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
print("===================================================================")
print("month from " + month_s + " to " + month_e + ", ntime = " + ntime)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
t     = new((/nyear,ntime,nlev,nlat,nlon/),float)
w     = new((/nyear,ntime,nlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("handle with " + " " + years(nt))
nv = 0
f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
printVarSummary(w)

if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb    
    opt = False
    edd_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    edd_v = bw_bandpass_filter(v,fca,fcb,opt,1)
else
    clm_u = conform(u, dim_avg_n(u,iyear), (/1,2,3,4/))
    clm_v = conform(v, dim_avg_n(v,iyear), (/1,2,3,4/)) 
    clm_t = conform(t, dim_avg_n(t,iyear), (/1,2,3,4/)) 
    clm_w = conform(w, dim_avg_n(w,iyear), (/1,2,3,4/)) 
    edd_u = u - clm_u
    edd_v = v - clm_v   
    edd_t = t - clm_t   
    edd_w = w - clm_w  
end if
printVarSummary(clm_u)
delete([/u,v,t,w/])
delete([/clm_u,clm_v,clm_t,clm_w/])

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
coslat = conform(edd_u,cos(lat),ilat)
Qeddy6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Qeddy6(:,0,:,:,:,:) = - (1.0/a/coslat) * center_finite_diff_n(edd_u*edd_t       ,lon,cyclic,opt,ilon)
Qeddy6(:,1,:,:,:,:) = - (1.0/a/coslat) * center_finite_diff_n(edd_v*edd_t*coslat,lat,cyclic,opt,ilat) 
Qeddy6(:,2,:,:,:,:) = - center_finite_diff_n(edd_w*edd_t,lev*100,cyclic,opt,ilev) 
Qeddy = tofloat(dim_avg_n(Qeddy6,iyear))
Qeddy_m(:,nm,:,:,:) = dim_avg_n(Qeddy,itime)
print(Qeddy(0,15,5,30,30:35))
delete([/Qeddy6,Qeddy/])

sinlat = conform(edd_u,sin(lat),ilat)
A6   = new((/nyear,2,ntime,nlev,nlat,nlon/),double)
vort = (center_finite_diff_n(edd_v,lon,cyclic,opt,ilon) - coslat*center_finite_diff_n(edd_u,lat,cyclic,opt,ilat) + sinlat*edd_u )/a/coslat
A6(:,0,:,:,:,:) = -(1.0/a/a/coslat/coslat)*center_finite_diff_n(edd_u*vort, lon, cyclic, opt, ilon)
A6(:,1,:,:,:,:) = -(1.0/a/a/coslat       )*center_finite_diff_n(edd_v*vort, lat, cyclic, opt, ilat)
A = tofloat(dim_avg_n(A6,iyear))
A_m(:,nm,:,:,:) = dim_avg_n(A,itime)
print(A(0,15,5,30,30:35))
print(vort(0,15,5,30,30:35))
delete([/A6,A,coslat,sinlat/])

Sflux6   = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Sflux6(:,0,:,:,:,:) = edd_u*vort
Sflux6(:,1,:,:,:,:) = edd_v*vort
Sflux6(:,2,:,:,:,:) = edd_w*vort
Sflux = tofloat(dim_avg_n(Sflux6,iyear))
Sflux_m(:,nm,:,:,:) = dim_avg_n(Sflux,itime)
print(Sflux(0,15,5,30,30:35))
delete([/Sflux6,Sflux,vort/])

Tflux6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Tflux6(:,0,:,:,:,:) = edd_u*edd_t
Tflux6(:,1,:,:,:,:) = edd_v*edd_t
Tflux6(:,2,:,:,:,:) = edd_w*edd_t
Tflux = tofloat(dim_avg_n(Tflux6,iyear))
Tflux_m(:,nm,:,:,:) = dim_avg_n(Tflux,itime)
print(Tflux(0,15,5,30,30:35))
delete([/Tflux6,Tflux/])

Mflux6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Mflux6(:,0,:,:,:,:) = edd_u*edd_v
Mflux6(:,1,:,:,:,:) = edd_u*edd_w
Mflux6(:,2,:,:,:,:) = edd_w*edd_v
Mflux = tofloat(dim_avg_n(Mflux6,iyear))
Mflux_m(:,nm,:,:,:) = dim_avg_n(Mflux,itime)
print(Mflux(0,15,5,30,30:35))
delete([/Mflux6,Mflux/])
delete([/edd_u,edd_v,edd_t,edd_w/])
end do 

Qeddy@long_name = "four terms of transient eddy heating, -dutdx, -dvtdy, -dwtdz, wt"
Qeddy@units     = "K/s" 
A@long_name = "two terms of diabatic heating, -duvortdx, -dvvortdy"
A@units     = "s^-2" 
Tflux@long_name = "three terms of transient eddy heating, ut, vt, wt"
Tflux@units     = "K*m/s" 
Sflux@long_name = "three terms of vorticity transport, uvort, vvort, wvort"
Sflux@units     = "m*s^-2" 
Mflux@long_name = "three terms of momentum transport, uv, uw, vw"
Mflux@units     = "m^2*s^-2" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "1-8 monthly three forcing factors, Qd, Qeddy, A"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"term1","term2","time","lev","lat","lon"/);
dimSizes = (/ 3     , 2     ,nmonth,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False  ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"Qeddy",typeof(Qeddy_m),(/"term1","time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A_m),(/"term2","time","lev","lat","lon"/))
filevardef(ncdf,"Tflux",typeof(Tflux_m),(/"term1","time","lev","lat","lon"/))
filevardef(ncdf,"Sflux",typeof(Sflux_m),(/"term1","time","lev","lat","lon"/))
filevardef(ncdf,"Mflux",typeof(Mflux_m),(/"term1","time","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"Qeddy",Qeddy_m)
filevarattdef(ncdf,"A"    ,A_m)
filevarattdef(ncdf,"Tflux",Tflux_m)
filevarattdef(ncdf,"Sflux",Sflux_m)
filevarattdef(ncdf,"Mflux",Mflux_m)

ncdf->time = (/12,34,56,78/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->Qeddy = (/Qeddy_m/)
ncdf->A     = (/A_m/)
ncdf->Tflux = (/Tflux_m/)
ncdf->Sflux = (/Sflux_m/)
ncdf->Mflux = (/Mflux_m/)
end

