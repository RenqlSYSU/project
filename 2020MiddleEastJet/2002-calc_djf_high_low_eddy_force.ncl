;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1  = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","omega","air","hgt"/)
fileout  = path2 + "ERA-Interim_DJF_Qeddy_A_freq-year.nc"

year_option = True;False 
geo = True ;False ;if true, only can calc 10N 

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

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
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)
nlev  = dimsizes(lev)
znlev = 5
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
w     := new((/nyear,ntime,nlev,nlat,nlon/),float)
z     := new((/nyear,ntime,znlev,nlat,nlon/),float)

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
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 4
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev(0:(znlev-1))},{lats:latn},{lonl:lonr})
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
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 4
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev(0:(znlev-1))},{lats:latn},{lonl:lonr})
end if
end do

if (any(ismissing(u))) then
    print("uwnd has missing")
end if
if (any(ismissing(v))) then
    print("vwnd has missing")
end if
if (any(ismissing(w))) then
    print("omega has missing")
end if
if (any(ismissing(t))) then
    print("air has missing")
end if
if (any(ismissing(z))) then
    print("hgt has missing")
end if
z = linmsg_n(z,-1,ilon)
;--------------------------------------
;The underlying wind field is replaced by quasi-geostrophic relations
;---------------------------------------------------
if(geo) then
dlev = -75*100 ;Pa
coslat = conform(z(:,:,0,:,:),cos(lat),ilat-1)
f0     = conform(z(:,:,0,:,:),2*(2*pi/24.0/3600.0)*sin(lat),ilat-1)
;z = z*g
do nlv = 0,1,1
    nl = 1 - nlv
    z1 = z(:,:,nl+2,:,:) + 2*R*t(:,:,nl+1,:,:)*dlev/lev(nl+1)/100
    u1 = -center_finite_diff_n(z(:,:,nl,:,:),lat,cyclic,opt,ilat-1)/f0/a
    v1 =  center_finite_diff_n(z(:,:,nl,:,:),lon,loncyc,opt,ilon-1)/f0/a/coslat
    t1 = -lev(nl)*100*(z(:,:,nl+1,:,:)-z(:,:,nl,:,:))/dlev/R
    print(max(abs(z(:,:,nl,:,:)-z1)))
    print(max(abs(u(:,:,nl,:,:)-u1)))
    print(max(abs(v(:,:,nl,:,:)-v1)))
    print(max(abs(t(:,:,nl,:,:)-t1)))
    print(z(1,1,nl,10,10:20)-z1(1,1,10,10:20))
    print(u(1,1,nl,10,10:20)-u1(1,1,10,10:20))
    print(v(1,1,nl,10,10:20)-v1(1,1,10,10:20))
    print(t(1,1,nl,10,10:20)-t1(1,1,10,10:20))
    ;z(:,:,nl,:,:) = (/z1/)
    u(:,:,nl,:,:) = (/u1/)
    v(:,:,nl,:,:) = (/v1/)
    ;t(:,:,nl,:,:) = (/t1/)
end do
delete([/coslat,f0,dlev,u1,v1,z1,t1,z/])
end if
printVarSummary(u)

if (any(ismissing(u))) then
    print("uwnd has missing")
end if
if (any(ismissing(v))) then
    print("vwnd has missing")
end if

;--------------------------------------
;DJF
;---------------------------------------------------
u2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
t2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
w2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    t2(ny,0:30,:,:,:) = t(ny,59:89,:,:,:)
    w2(ny,0:30,:,:,:) = w(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    t2(ny,31:89,:,:,:) = t(ny+1,0:58,:,:,:)
    w2(ny,31:89,:,:,:) = w(ny+1,0:58,:,:,:)
end do
u := u2
v := v2
t := t2
w := w2
printVarSummary(u)

;--------------------------------------
;The data of 15 layers is inserted into 21 layers by log interpolation
;---------------------------------------------------
lev0 = lev
lev := (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,175,150,125, 100/)
nlev = dimsizes(lev)
Qeddh := new((/nyear-1,ntime,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(u    ,lev0,ilev)
p_out  = conform(Qeddh,lev ,ilev)
u := int2p_n(p_in,u,p_out,linlog,ilev)
v := int2p_n(p_in,v,p_out,linlog,ilev)
t := int2p_n(p_in,t,p_out,linlog,ilev)
w := int2p_n(p_in,w,p_out,linlog,ilev)
delete([/p_in,p_out,linlog/])

;--------------------------------------
;Calculate seasonal averages and anomalies based on seasonal averages
;---------------------------------------------------
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

cf = ezfftf_n(t,itime)
term = cf 
term(:,:,hgh_cf,:,:,:) = 0.0  
low_t = ezfftb_n(term,0,itime+1)
term = cf 
term(:,:,low_cf,:,:,:) = 0.0  
hgh_t = ezfftb_n(term,0,itime+1)
printVarSummary(hgh_t)

cf = ezfftf_n(w,itime)
term = cf 
term(:,:,hgh_cf,:,:,:) = 0.0  
low_w = ezfftb_n(term,0,itime+1)
term = cf 
term(:,:,low_cf,:,:,:) = 0.0  
hgh_w = ezfftb_n(term,0,itime+1)
printVarSummary(hgh_w)
delete([/cf,term,fcf,fca,hgh_cf,low_cf/])

clm_t = dim_avg_n(t,itime) 
delete([/u,v,t,w/])

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
lev4d = conform(clm_t, lev*100,ilev-1)
sigm0 = (R/cp)*(clm_t/lev4d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev-1)
delete([/lev4d,clm_t/])

coslat = conform(Qeddh,cos(lat),ilat)
tanlat = conform(Qeddh,tan(lat),ilat)
lev5d  = conform(Qeddh, lev*100,ilev)

Qeddh := - (1.0/a/coslat) * center_finite_diff_n(hgh_u*hgh_t       ,lon,loncyc,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(hgh_v*hgh_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(hgh_w*hgh_t)/lev5d - center_finite_diff_n(hgh_w*hgh_t,lev*100,cyclic,opt,ilev) 
Qeddh := dim_avg_n(Qeddh,itime)
printVarSummary(Qeddh)
print(Qeddh(0,5,30,30))

Qeddl := - (1.0/a/coslat) * center_finite_diff_n(low_u*low_t       ,lon,loncyc,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(low_v*low_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(low_w*low_t)/lev5d - center_finite_diff_n(low_w*low_t,lev*100,cyclic,opt,ilev) 
Qeddl := dim_avg_n(Qeddl,itime)
printVarSummary(Qeddl)
print(Qeddl(0,5,30,30))
delete([/lev5d/])

Ah:= -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_v              ,lon    ,loncyc,opt,ilon),lon    ,loncyc,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_v*hgh_v*coslat       ,lon    ,loncyc,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(hgh_v*hgh_w              ,lon    ,loncyc,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(hgh_u*hgh_u              ,lon    ,loncyc,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_u              ,lat    ,cyclic,opt,ilat),lon    ,loncyc,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(hgh_u*hgh_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     
     -(tanlat/a                )*                     center_finite_diff_n(hgh_u*hgh_w              ,lev*100,cyclic,opt,ilev)
Ah := dim_avg_n(Ah,itime)
printVarSummary(Ah)
print(Ah(0,5,30,30))
delete([/hgh_u,hgh_v,hgh_t,hgh_w/])

Al:= -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(low_u*low_v              ,lon    ,loncyc,opt,ilon),lon    ,loncyc,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(low_v*low_v*coslat       ,lon    ,loncyc,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(low_v*low_w              ,lon    ,loncyc,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(low_u*low_u              ,lon    ,loncyc,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(low_u*low_u              ,lat    ,cyclic,opt,ilat),lon    ,loncyc,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(low_u*low_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(low_u*low_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     +(1.0   /a                )*center_finite_diff_n(center_finite_diff_n(low_u*low_w              ,lat    ,cyclic,opt,ilat),lev*100,cyclic,opt,ilev)\
     -(tanlat/a                )*                     center_finite_diff_n(low_u*low_w              ,lev*100,cyclic,opt,ilev)
Al := dim_avg_n(Al,itime)
printVarSummary(Al)
print(Al(0,5,30,30))
delete([/low_u,low_v,low_t,low_w/])
delete([/coslat,tanlat/])

Qeddh@long_name = "Climatology 2-10-day transient eddy heating"
Qeddh@units     = "K/s" 
Qeddl@long_name = "Climatology 10-90-day transient eddy heating"
Qeddl@units     = "K/s" 
Ah@long_name = "Climatology 2-10-day diabatic heating"
Ah@units     = "s^-2" 
Al@long_name = "Climatology 10-90-day diabatic heating"
Al@units     = "s^-2" 

if(.not.(year_option)) then 
    sigm0 := dim_avg_n_Wrap(sigm0,0)
    Qeddh := dim_avg_n_Wrap(Qeddh,0)
    Ah    := dim_avg_n_Wrap(Ah,0)
    Qeddl := dim_avg_n_Wrap(Qeddl,0)
    Al    := dim_avg_n_Wrap(Al,0)
end if

;==========================================================
;calc Qd foring and Qeddy forcing, dimension is 4d 
;============================================================
if(year_option) then 
lev4d = conform(Ah, lev*100,ilev-1)
sigm1 = sigm0*R/lev4d
sigm1_m = conform(lev4d, wgt_areaave(sigm1,1.0,1.0,0), (/iyear,ilev-1/))
f0    = conform(Ah, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-1)
f_Qeddh  = -f0*R*center_finite_diff_n((Qeddh/sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
f_Qeddl  = -f0*R*center_finite_diff_n((Qeddl/sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
delete(lev4d)
else
lev3d = conform(Ah, lev*100,ilev-2)
sigm1 = sigm0*R/lev3d
sigm1_m = conform(lev3d, wgt_areaave(sigm1,1.0,1.0,0), (/ilev-2/))
f0    = conform(Ah, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-2)
f_Qeddh  = -f0*R*center_finite_diff_n((Qeddh/sigm1_m/lev3d),lev*100,cyclic,opt,ilev-2)
f_Qeddl  = -f0*R*center_finite_diff_n((Qeddl/sigm1_m/lev3d),lev*100,cyclic,opt,ilev-2)
delete(lev3d)
end if
delete(f0)

f_Qeddh@long_name = "Climatology 2-10-day transient eddy heating"
f_Qeddh@units     = "s^-2"                 
f_Qeddl@long_name = "Climatology 10-90-day transient eddy heating"
f_Qeddl@units     = "s^-2" 
sigm1@long_name = "One kind of static stability measured by potential temperature"
sigm1@units     = "J/(Pa^2*Kg)" 
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
fileatt@discription   = "DJF mean four forces, Qeddh, Ah ,Qeddl, Al"
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
if(year_option) then 
filevardef(ncdf,"sigm1",typeof(sigm1),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Qeddl",typeof(Qeddl),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Al"    ,typeof(Al),(/"year","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddl",typeof(f_Qeddl),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Qeddh",typeof(Qeddh),(/"year","lev","lat","lon"/))
filevardef(ncdf,"Ah"    ,typeof(Ah),(/"year","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddh",typeof(f_Qeddh),(/"year","lev","lat","lon"/))
else
filevardef(ncdf,"sigm1",typeof(sigm1),(/"lev","lat","lon"/))
filevardef(ncdf,"Qeddl",typeof(Qeddl),(/"lev","lat","lon"/))
filevardef(ncdf,"Al"    ,typeof(Al),(/"lev","lat","lon"/))
filevardef(ncdf,"f_Qeddl",typeof(f_Qeddl),(/"lev","lat","lon"/))
filevardef(ncdf,"Qeddh",typeof(Qeddh),(/"lev","lat","lon"/))
filevardef(ncdf,"Ah"    ,typeof(Ah),(/"lev","lat","lon"/))
filevardef(ncdf,"f_Qeddh",typeof(f_Qeddh),(/"lev","lat","lon"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm1",sigm1)
filevarattdef(ncdf,"Qeddl",Qeddl)
filevarattdef(ncdf,"Al"    ,Al)
filevarattdef(ncdf,"f_Qeddl",f_Qeddl)
filevarattdef(ncdf,"Qeddh",Qeddh)
filevarattdef(ncdf,"Ah"    ,Ah)
filevarattdef(ncdf,"f_Qeddh",f_Qeddh)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm1 = (/sigm1/)
ncdf->Qeddl = (/Qeddl/)
ncdf->Al     = (/Al/)
ncdf->f_Qeddl = (/f_Qeddl/)
ncdf->Qeddh = (/Qeddh/)
ncdf->Ah     = (/Ah/)
ncdf->f_Qeddh = (/f_Qeddh/)
end

