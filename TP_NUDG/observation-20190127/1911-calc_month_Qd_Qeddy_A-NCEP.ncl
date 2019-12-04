;***************************************************************************
;use the post-processing data to calculate climatology daily transient eddies
;1. do bandpass filtering;  2. calc transient eddies kinetic energy (EkE)
;3. calc climatology daily EKE
;4. store EKE on each level
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/01/25
;******************************************************************************
begin
;path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","omega","air","hgt"/)
option   = 0 ;0 means remove seasonal mean, 1 means bandpass filter
geo = False
csa = False 

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

if(option.eq.0) then
    fileout  = path2 + "mdata/NCEP1_DJF_Qd_Qeddy_A_1204.nc"
else                                     
    fileout  = path2 + "mdata/ERA-Interim_DJF_Qd_Qeddy_A-" + cb + "-" + ca + "day.nc"
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
lats =  -2 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100/) 
;lev  = (/1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+"."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
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
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
select_time := ind(time(:,1).eq.12.or.time(:,1).le.2)
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
w     := new((/nyear,ntime,nlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    ;if(nm.gt.0) then 
    ;select_time = select_time + 1
    ;end if
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    ;nv = 4
    ;f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    ;z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    ;nv = 4
    ;f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    ;z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end if
end do
printVarSummary(u)

if(geo) then
dlev = -75*100 ;Pa
coslat = conform(z(:,:,0,:,:),cos(lat),ilat-1)
f0     = conform(z(:,:,0,:,:),2*(2*pi/24.0/3600.0)*sin(lat),ilat-1)
z = z*g
do nlv = 0,0,1
    nl = 0 - nlv
    z1 = z(:,:,nl+2,:,:) + 2*R*t(:,:,nl+1,:,:)*dlev/lev(nl+1)/100
    u1 = -center_finite_diff_n(z1,lat,cyclic,opt,ilat-1)/f0/a
    v1 =  center_finite_diff_n(z1,lon,loncyc,opt,ilon-1)/f0/a/coslat
    t1 = -lev(nl)*100*(z(:,:,nl+1,:,:)-z1)/dlev/R
    print(max(abs(z(:,:,nl,:,:)-z1)))
    print(max(abs(u(:,:,nl,:,:)-u1)))
    print(max(abs(v(:,:,nl,:,:)-v1)))
    print(max(abs(t(:,:,nl,:,:)-t1)))
    print(z(1,1,nl,10,10:20)-z1(1,1,10,10:20))
    print(u(1,1,nl,10,10:20)-u1(1,1,10,10:20))
    print(v(1,1,nl,10,10:20)-v1(1,1,10,10:20))
    print(t(1,1,nl,10,10:20)-t1(1,1,10,10:20))
    z(:,:,nl,:,:) = (/z1/)
    u(:,:,nl,:,:) = (/u1/)
    v(:,:,nl,:,:) = (/v1/)
    t(:,:,nl,:,:) = (/t1/)
end do
delete([/coslat,f0,dlev,u1,v1,z1,t1,z/])
end if

if(csa) then 
u!0 = "year"
v!0 = "year"
t!0 = "year"
w!0 = "year"
lev0 = lev
lev := fspan(1000,100,19)
wts = tofloat(lev0) ;weights
wts = 1.0
;wts(0) = 0.1 ;1000
;wts(1) = 0.2 ;925
;wts(2) = 0.5 ;850
knots = nlev - 2
smth = 0.0
nderiv = 0
u_o = csa1x(lev0,u(year|:,time|:,lat|:,lon|:,level|:),wts,knots,smth,nderiv,lev)
v_o = csa1x(lev0,v(year|:,time|:,lat|:,lon|:,level|:),wts,knots,smth,nderiv,lev)
t_o = csa1x(lev0,t(year|:,time|:,lat|:,lon|:,level|:),wts,knots,smth,nderiv,lev)
w_o = csa1x(lev0,w(year|:,time|:,lat|:,lon|:,level|:),wts,knots,smth,nderiv,lev)
printVarSummary(u_o)
u_o!0 = "year"
u_o!1 = "time"
u_o!2 = "lat"
u_o!3 = "lon"
u_o!4 = "lev"
u := u_o(year|:,time|:,lev|:,lat|:,lon|:)
u@_FillValue = default_fillvalue("float")
v_o!0 = "year"
v_o!1 = "time"
v_o!2 = "lat"
v_o!3 = "lon"
v_o!4 = "lev"
v := v_o(year|:,time|:,lev|:,lat|:,lon|:)
v@_FillValue = default_fillvalue("float")
t_o!0 = "year"
t_o!1 = "time"
t_o!2 = "lat"
t_o!3 = "lon"
t_o!4 = "lev"
t := t_o(year|:,time|:,lev|:,lat|:,lon|:)
t@_FillValue = default_fillvalue("float")
w_o!0 = "year"
w_o!1 = "time"
w_o!2 = "lat"
w_o!3 = "lon"
w_o!4 = "lev"
w := w_o(year|:,time|:,lev|:,lat|:,lon|:)
w@_FillValue = default_fillvalue("float")
nlev = dimsizes(lev)
Qeddy := new((/nyear,ntime,nlev,nlat,nlon/),double)
else

lev0 = lev
lev := fspan(1000,100,19)
nlev = dimsizes(lev)
Qeddy := new((/nyear,ntime,nlev,nlat,nlon/),double)
Qd    := new((/nyear,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(u    ,lev0,ilev)
p_out  = conform(Qeddy,lev ,ilev)
u := int2p_n(p_in,u,p_out,linlog,ilev)
v := int2p_n(p_in,v,p_out,linlog,ilev)
t := int2p_n(p_in,t,p_out,linlog,ilev)
w := int2p_n(p_in,w,p_out,linlog,ilev)
delete([/p_in,p_out,linlog/])
end if

clm_u = dim_avg_n(u,itime)
clm_v = dim_avg_n(v,itime) 
clm_t = dim_avg_n(t,itime) 
clm_w = dim_avg_n(w,itime) 
if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb    
    opt = False
    edd_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    edd_v = bw_bandpass_filter(v,fca,fcb,opt,1)
else
    edd_u = u - conform(u, clm_u, (/0,2,3,4/))
    edd_v = v - conform(v, clm_v, (/0,2,3,4/))   
    edd_t = t - conform(t, clm_t, (/0,2,3,4/))   
    edd_w = w - conform(w, clm_w, (/0,2,3,4/))  
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
print(Qd  (1,5,30,30))
delete([/clm_u,clm_v,clm_t,clm_w/])
delete([/coslat4d,tanlat4d,lev4d/])

coslat = conform(Qeddy,cos(lat),ilat)
tanlat = conform(Qeddy,tan(lat),ilat)
lev5d  = conform(Qeddy, lev*100,ilev)

sig0e = (R/cp)*(edd_t/lev5d) - center_finite_diff_n(edd_t,lev*100,cyclic,opt,ilev)
Qeddy = - (edd_u/a/coslat) * center_finite_diff_n(edd_t,lon,loncyc,opt,ilon) \
        - (edd_v/a       ) * center_finite_diff_n(edd_t,lat,cyclic,opt,ilat) \
        + sig0e * edd_w 
;Qeddy := - (1.0/a/coslat) * center_finite_diff_n(edd_u*edd_t       ,lon,loncyc,opt,ilon) \
;         - (1.0/a/coslat) * center_finite_diff_n(edd_v*edd_t*coslat,lat,cyclic,opt,ilat) \
;         + (R/cp)*(edd_w*edd_t)/lev5d - center_finite_diff_n(edd_w*edd_t,lev*100,cyclic,opt,ilev) 

Qeddy := dim_avg_n(Qeddy,itime)
print(Qeddy(1,5,30,30))
Qd   = Qd - Qeddy
delete([/sig0e,lev5d/])

A := -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v              ,lon    ,loncyc,opt,ilon),lon    ,loncyc,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_v*edd_v*coslat       ,lon    ,loncyc,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(edd_v*edd_w              ,lon    ,loncyc,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(edd_u*edd_u              ,lon    ,loncyc,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_u              ,lat    ,cyclic,opt,ilat),lon    ,loncyc,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     +(1.0   /a                )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_w              ,lat    ,cyclic,opt,ilat),lev*100,cyclic,opt,ilev)\
     -(tanlat/a                )*                     center_finite_diff_n(edd_u*edd_w              ,lev*100,cyclic,opt,ilev)
A := dim_avg_n(A,itime)
print(A(1,5,30,30))
delete([/edd_u,edd_v,edd_t,edd_w/])
delete([/coslat,tanlat/])

Qeddy@long_name = "Climatology daily transient eddy heating"
Qeddy@units     = "K/s" 
Qd@long_name = "Climatology daily diabatic heating"
Qd@units     = "K/s" 
A@long_name = "Climatology daily diabatic heating"
A@units     = "s^-2" 

;==========================================================
;calc Qd foring and Qeddy forcing, dimension is 4d 
;============================================================
lev4d = conform(A, lev*100,ilev-1)
sigm1 = sigm0*R/lev4d
sigm1_m = conform(lev4d, wgt_areaave(sigm1,1.0,1.0,0), (/iyear,ilev-1/))
f0    = conform(A, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-1)

;do ilat = 0, nlat-1
;if (abs(lat(ilat)*180.0/pi ).lt. 10. ) then
;f0(:,:,ilat,:)= A@_FillValue
;end if
;end do

f_Qd     = -f0*R*center_finite_diff_n((Qd   /sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
f_Qeddy  = -f0*R*center_finite_diff_n((Qeddy/sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
delete(lev4d)
delete(f0)

f_Qeddy@long_name = "Climatology daily transient eddy heating"
f_Qeddy@units     = "s^-2" 
f_Qd@long_name    = "Climatology daily diabatic heatinig"
f_Qd@units        = "s^-2" 
sigm1@long_name = "One kind of static stability measured by potential temperature"
sigm1@units     = "J/(Pa^2*Kg)" 
print(f_Qeddy(1,5,30,30))
print(f_Qd   (1,5,30,30))
print(lat(0:9))

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
fileatt@discription   = "DJF mean three forcing factors, Qd, Qeddy, A"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year","lev","lat","lon"/);
dimSizes = (/nyear ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"sigm1",typeof(sigm1),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Qeddy",typeof(Qeddy),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Qd"   ,typeof(Qd),(/"year","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A),(/"year","lev","lat","lon"/))
filevardef(ncdf,"f_Qd",typeof(f_Qd),(/"year","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddy",typeof(f_Qeddy),(/"year","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm1",sigm1)
filevarattdef(ncdf,"Qeddy",Qeddy)
filevarattdef(ncdf,"Qd"   ,Qd)
filevarattdef(ncdf,"A"    ,A)
filevarattdef(ncdf,"f_Qeddy",f_Qeddy)
filevarattdef(ncdf,"f_Qd"   ,f_Qd)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm1 = (/sigm1/)
ncdf->Qeddy = (/Qeddy/)
ncdf->Qd    = (/Qd/)
ncdf->A     = (/A/)
ncdf->f_Qeddy = (/f_Qeddy/)
ncdf->f_Qd    = (/f_Qd/)
end

