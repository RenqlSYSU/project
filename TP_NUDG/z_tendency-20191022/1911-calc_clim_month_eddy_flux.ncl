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
path    = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
nmonth = 4
month_s = 1 
month_e = 8

;parameter for bw_bandpass_filter
ca = 7.0 ;day
cb = 2.5 

if(case_option.eq.0) then
    case    = (/"CTRL","NG24h","NG6h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_eddy_flux4.nc"
    else
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_eddy_flux_" + cb + "-" + ca + "day.nc"
    end if  
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux.nc"
    else                                                                                              
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux_" + cb + "-" + ca + "day.nc"
    end if  
end if

if(case_option.eq.2) then
    case    = (/"F2000","NUDG","NG58"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux12.nc"
    else
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux12_" + cb + "-" + ca + "day.nc"
    end if  
end if

var_name= (/"U","V","T","OMEGA"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats = 0   ;-90 ;
latn = 90  ; 90 ;
lonl = 0   ; 0  ;
lonr = 360 ;360 ;
;lev  = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels, hPa
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;13 levels ,hPa

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
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

do nf =2,2,1
;do nf = 0,ncase-2,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nm = 0, nmonth-1, 1
month1_s = nm*2 + 1
month1_e = nm*2 + 2
select_time := ind(time(:,1).ge.month1_s.and.time(:,1).le.month1_e)
ntime = dimsizes(select_time)
print("===================================================================")
print("month from " + month1_s + " to " + month1_e + ", ntime = " + ntime)
u    = new((/nyear,ntime,nlev,nlat,nlon/),float)
v    = new((/nyear,ntime,nlev,nlat,nlon/),float)
t    = new((/nyear,ntime,nlev,nlat,nlon/),float)
w    = new((/nyear,ntime,nlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("handle with " + case(nf) + " " + years(nt))
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
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
lev5d = tofloat(conform(clm_t,lev*100,ilev))
density = 1.293*(lev5d/101325.0)*(273.15/clm_t)
delete([/u,v,t,w/])
delete([/clm_u,clm_v,clm_t,clm_w,lev5d/])

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
A6(:,0,:,:,:,:) = -(1.0/a/coslat)*center_finite_diff_n(edd_u*vort       , lon, cyclic, opt, ilon)
A6(:,1,:,:,:,:) = -(1.0/a/coslat)*center_finite_diff_n(edd_v*vort*coslat, lat, cyclic, opt, ilat)
A = tofloat(dim_avg_n(A6,iyear))
A_m(:,nm,:,:,:) = dim_avg_n(A,itime)
print(A(0,15,5,30,30:35))
print(vort(0,15,5,30,30:35))
delete([/A6,A,coslat,sinlat/])

Sflux6   = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Sflux6(:,0,:,:,:,:) = edd_u*vort
Sflux6(:,1,:,:,:,:) = edd_v*vort
Sflux6(:,2,:,:,:,:) = edd_w*vort/density/g
Sflux = tofloat(dim_avg_n(Sflux6,iyear))
Sflux_m(:,nm,:,:,:) = dim_avg_n(Sflux,itime)
print(Sflux(0,15,5,30,30:35))
delete([/Sflux6,Sflux,vort/])

Tflux6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Tflux6(:,0,:,:,:,:) = edd_u*edd_t
Tflux6(:,1,:,:,:,:) = edd_v*edd_t
Tflux6(:,2,:,:,:,:) = edd_w*edd_t/density/g
Tflux = tofloat(dim_avg_n(Tflux6,iyear))
Tflux_m(:,nm,:,:,:) = dim_avg_n(Tflux,itime)
print(Tflux(0,15,5,30,30:35))
delete([/Tflux6,Tflux/])

Mflux6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Mflux6(:,0,:,:,:,:) = edd_u*edd_v
Mflux6(:,1,:,:,:,:) = edd_u*edd_w/density/g
Mflux6(:,2,:,:,:,:) = edd_w*edd_v/density/g
Mflux = tofloat(dim_avg_n(Mflux6,iyear))
Mflux_m(:,nm,:,:,:) = dim_avg_n(Mflux,itime)
print(Mflux(0,15,5,30,30:35))
delete([/Mflux6,Mflux/])
delete([/edd_u,edd_v,edd_t,edd_w,density/])
end do 
Qeddy_m@long_name = "four terms of transient eddy heating, -dutdx, -dvtdy, -dwtdz, wt"
Qeddy_m@units     = "K/s" 
A_m@long_name = "two terms of diabatic heating, -duvortdx, -dvvortdy"
A_m@units     = "s^-2" 
Tflux_m@long_name = "three terms of transient eddy heating, ut, vt, wt"
Tflux_m@units     = "K*m/s" 
Sflux_m@long_name = "three terms of vorticity transport, uvort, vvort, wvort"
Sflux_m@units     = "m*s^-2" 
Mflux_m@long_name = "three terms of momentum transport, uv, uw, vw"
Mflux_m@units     = "m^2*s^-2" 

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
fileatt@discription   = "12 month"
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

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->Qeddy = (/Qeddy_m/)
ncdf->A     = (/A_m/)
ncdf->Tflux = (/Tflux_m/)
ncdf->Sflux = (/Sflux_m/)
ncdf->Mflux = (/Mflux_m/)
end do
end

