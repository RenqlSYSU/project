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
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_Qd_Qeddy_A-g_mis.nc"
    else
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_daily_Qd_Qeddy_A_" + cb + "-" + ca + "day.nc"
    end if  
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_Qd_Qeddy_A.nc"
    else                                                      
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_Qd_Qeddy_A_" + cb + "-" + ca + "day.nc"
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
lats = -2 ;-90 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time = cd_calendar(f->time,0)
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
Qd_m       = new((/4,nlev,nlat,nlon/),double)
Qeddy_m    = new((/4,nlev,nlat,nlon/),double)
A_m        = new((/4,nlev,nlat,nlon/),double)
f_Qd_m     = new((/4,nlev,nlat,nlon/),double)
f_Qeddy_m  = new((/4,nlev,nlat,nlon/),double)
sigm0_m    = new((/4,nlev,nlat,nlon/),double)

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
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 0,2,2
if((case_option.eq.1).and.(nf.eq.1))
    years  = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nm = 0, 3, 1
month_s = 2*nm + 1
month_e = 2*nm + 2
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
w     := new((/nyear,ntime,nlev,nlat,nlon/),float)
Qeddy := new((/nyear,ntime,nlev,nlat,nlon/),double)

do nt = 0,nyear-1,1
print("handle with " + case(nf) + " " + years(nt))
nv = 0
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
;u = linmsg_n(u,-1,ilon)
;v = linmsg_n(v,-1,ilon)
;t = linmsg_n(t,-1,ilon)
;w = linmsg_n(w,-1,ilon)
printVarSummary(w)

if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb    
    opt = False
    edd_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    edd_v = bw_bandpass_filter(v,fca,fcb,opt,1)
else
    clm_u = dim_avg_n(u,iyear)
    clm_v = dim_avg_n(v,iyear) 
    clm_t = dim_avg_n(t,iyear) 
    clm_w = dim_avg_n(w,iyear) 
    edd_u = u - conform(u, clm_u, (/1,2,3,4/))
    edd_v = v - conform(v, clm_v, (/1,2,3,4/))   
    edd_t = t - conform(t, clm_t, (/1,2,3,4/))   
    edd_w = w - conform(w, clm_w, (/1,2,3,4/))  
end if
printVarSummary(clm_u)
delete([/u,v,t,w/])

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
coslat4d = conform(clm_u,cos(lat),ilat-1)
tanlat4d = conform(clm_u,tan(lat),ilat-1)
lev4d    = conform(clm_u, lev*100,ilev-1)

sigm0 = (R/cp)*(clm_t/lev4d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev-1)
Qd    =   (clm_u/a/coslat4d) * center_finite_diff_n(clm_t,lon,loncyc,opt,ilon-1) \
        + (clm_v/a         ) * center_finite_diff_n(clm_t,lat,cyclic,opt,ilat-1) \
        - sigm0 * clm_w 
delete([/clm_u,clm_v,clm_t,clm_w/])
delete([/coslat4d,tanlat4d,lev4d/])

coslat = conform(Qeddy,cos(lat),ilat)
tanlat = conform(Qeddy,tan(lat),ilat)
lev5d  = conform(Qeddy, lev*100,ilev)

Qeddy := - (1.0/a/coslat) * center_finite_diff_n(edd_u*edd_t       ,lon,loncyc,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(edd_v*edd_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(edd_w*edd_t)/lev5d - center_finite_diff_n(edd_w*edd_t,lev*100,cyclic,opt,ilev) 
Qeddy := dim_avg_n(Qeddy,iyear)
Qd   = Qd - Qeddy
print(Qeddy(15,5,30,30))
print(Qd  (15,5,30,30))
Qeddy_m(nm,:,:,:) = dim_avg_n(Qeddy,itime-1)
Qd_m(nm,:,:,:)    = dim_avg_n(Qd,itime-1)
delete(Qeddy)
delete(Qd)

A := -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v              ,lon    ,loncyc,opt,ilon),lon    ,loncyc,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_v*edd_v*coslat       ,lon    ,loncyc,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(edd_v*edd_w              ,lon    ,loncyc,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(edd_u*edd_u              ,lon    ,loncyc,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_u              ,lat    ,cyclic,opt,ilat),lon    ,loncyc,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     +(1.0   /a                )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_w              ,lat    ,cyclic,opt,ilat),lev*100,cyclic,opt,ilev)\
     -(tanlat/a                )*                     center_finite_diff_n(edd_u*edd_w              ,lev*100,cyclic,opt,ilev)
A := dim_avg_n(A,iyear)
print(A(15,5,30,30))
A_m(nm,:,:,:) = dim_avg_n(A,itime-1)
sigm0_m(nm,:,:,:) = dim_avg_n(sigm0,itime-1)
delete([/edd_u,edd_v,edd_t,edd_w/])
delete([/coslat,tanlat,lev5d,A,sigm0/])
end do

Qeddy_m@long_name = "Climatology daily transient eddy heating"
Qeddy_m@units     = "K/s" 
Qd_m@long_name = "Climatology daily diabatic heating"
Qd_m@units     = "K/s" 
A_m@long_name = "Climatology daily diabatic heating"
A_m@units     = "s^-2" 

;==========================================================
;calc Qd foring and Qeddy forcing, dimension is 4d 
;============================================================
lev4d = conform(A_m, lev*100,ilev-1)
sigm1 = sigm0_m*R/lev4d
f0    = conform(A_m, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-1)

;do ilat = 0, nlat-1
;if (abs(lat(ilat)*180.0/pi ).lt. 10. ) then
;f0(:,:,ilat,:)= A@_FillValue
;end if
;end do

f_Qd_m     = -f0*R*center_finite_diff_n((Qd_m   /sigm1/lev4d),lev*100,cyclic,opt,ilev-1)
f_Qeddy_m  = -f0*R*center_finite_diff_n((Qeddy_m/sigm1/lev4d),lev*100,cyclic,opt,ilev-1)
delete(lev4d)
delete(f0)

f_Qeddy_m@long_name = "Climatology daily transient eddy heating"
f_Qeddy_m@units     = "s^-2" 
f_Qd_m@long_name    = "Climatology daily diabatic heatinig"
f_Qd_m@units        = "s^-2" 
sigm1@long_name = "One kind of static stability measured by potential temperature"
sigm1@units     = "J/(Pa^2*Kg)" 
print(f_Qeddy_m(1,5,30,30))
print(f_Qd_m   (1,5,30,30))
print(lat(0:9))

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
dimNames = (/"time","lev","lat","lon"/);
dimSizes = (/   4  ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"sigm1",typeof(sigm1),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qeddy",typeof(Qeddy_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qd"   ,typeof(Qd_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qd",typeof(f_Qd_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddy",typeof(f_Qeddy_m),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm1",sigm1)
filevarattdef(ncdf,"Qeddy",Qeddy_m)
filevarattdef(ncdf,"Qd"   ,Qd_m)
filevarattdef(ncdf,"A"    ,A_m)
filevarattdef(ncdf,"f_Qeddy",f_Qeddy_m)
filevarattdef(ncdf,"f_Qd"   ,f_Qd_m)

ncdf->time = (/12,34,56,78/) 
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm1 = (/sigm1/)
ncdf->Qeddy = (/Qeddy_m/)
ncdf->Qd    = (/Qd_m/)
ncdf->A     = (/A_m/)
ncdf->f_Qeddy = (/f_Qeddy_m/)
ncdf->f_Qd    = (/f_Qd_m/)
end do
end

