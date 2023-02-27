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
var_name = (/"uwnd","vwnd","omega","air","hgt"/)

year_option = True;False ;
geo = True ;False ;if true, only can calc 10N 
fileout  = path2 + "ERA-Interim_DJF_in3-year.nc"

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
ilev  = 1
ilat  = 2
ilon  = 3
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
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 4
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev(0:(znlev-1))},{lats:latn},{lonl:lonr})
end if
end do

;--------------------------------------
;The underlying wind field is replaced by quasi-geostrophic relations
;---------------------------------------------------
if(geo) then
dlev = -75*100 ;Pa
coslat = conform(z(:,:,0,:,:),cos(lat),ilat)
f0     = conform(z(:,:,0,:,:),2*(2*pi/24.0/3600.0)*sin(lat),ilat)
;z = z*g
do nlv = 0,1,1
    nl = 1 - nlv
    z1 = z(:,:,nl+2,:,:) + 2*R*t(:,:,nl+1,:,:)*dlev/lev(nl+1)/100
    u1 = -center_finite_diff_n(z(:,:,nl,:,:),lat,cyclic,opt,ilat)/f0/a
    v1 =  center_finite_diff_n(z(:,:,nl,:,:),lon,loncyc,opt,ilon)/f0/a/coslat
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

;--------------------------------------
;DJF
;---------------------------------------------------
u2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
t2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    t2(ny,0:30,:,:,:) = t(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    t2(ny,31:89,:,:,:) = t(ny+1,0:58,:,:,:)
end do
u := u2
v := v2
t := t2
printVarSummary(u)

clm_u = dim_avg_n(u,1)
clm_v = dim_avg_n(v,1) 
clm_t = dim_avg_n(t,1) 
printVarSummary(clm_u)
delete([/u,v,t/])

;--------------------------------------
;The data of 15 layers is inserted into 21 layers by log interpolation
;---------------------------------------------------
lev0 = lev
lev := (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,175,150,125, 100/)
nlev = dimsizes(lev)
b1  := new((/nyear-1,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(clm_u,lev0,ilev)
p_out  = conform(b1   ,lev ,ilev)
clm_u := int2p_n(p_in,clm_u,p_out,linlog,ilev)
clm_v := int2p_n(p_in,clm_v,p_out,linlog,ilev)
clm_t := int2p_n(p_in,clm_t,p_out,linlog,ilev)
delete([/p_in,p_out,linlog/])

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
coslat = conform(clm_u,cos(lat),ilat)
lev4d  = conform(clm_u, lev*100,ilev)
sigm0  = (R/cp)*(clm_t/lev4d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev)
delete(lev4d)

;vort = (1/a/coslat)*(center_finite_diff_n(clm_v,lon,loncyc,opt,ilon) - center_finite_diff_n(clm_u*coslat,lat,cyclic,opt,ilat))
vort = uv2vr_cfd(clm_u,clm_v,vars&lat,vars&lon,3) 
b1x = - (clm_u/a/coslat)*center_finite_diff_n(vort,lon,loncyc,opt,ilon) 
b1y = - (clm_v/a       )*center_finite_diff_n(vort,lat,cyclic,opt,ilat)
print(vort(0,9,10,20:30))
print(b1x(0,9,10,20:30))

sigm0 = conform(clm_u,wgt_areaave(sigm0,1.0,1.0,0),(/iyear,ilev/))
f0   = conform(clm_u,2*(2*pi/24.0/3600.0)*sin(lat),ilat)

b2 = - 2*(2.0*pi/24.0/3600.0)*coslat*clm_v/a

t3x = -(clm_u/a/coslat)*center_finite_diff_n(clm_t,lon,loncyc,opt,ilon)
t3y = -(clm_v/a       )*center_finite_diff_n(clm_t,lat,cyclic,opt,ilat)
b3x = -f0 * center_finite_diff_n(t3x/sigm0,lev*100,cyclic,opt,ilev)
b3y = -f0 * center_finite_diff_n(t3y/sigm0,lev*100,cyclic,opt,ilev)
delete([/clm_u,clm_v,clm_t,coslat/])
delete([/f0/])

b2@long_name = "Horizontal advection of absolute vorticity"
b2@units     = "s^-2" 
b1x@long_name = "x Horizontal advection of relative vorticity"
b1x@units     = "s^-2" 
b3x@long_name = "Vertical gradient of x temperature horizontal advection"
b3x@units     = "s^-2" 
t3x@long_name = "x temperature horizontal advection"
t3x@units     = "K/s" 
b1y@long_name = "y Horizontal advection of relative vorticity"
b1y@units     = "s^-2" 
b3y@long_name = "Vertical gradient of y temperature horizontal advection"
b3y@units     = "s^-2" 
t3y@long_name = "y temperature horizontal advection"
t3y@units     = "K/s"

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
fileatt@discription   = "DJF daily three advection terms"
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
filevardef(ncdf,"b2",typeof(b2),(/"year","lev","lat","lon"/))
filevardef(ncdf,"b1x",typeof(b1x),(/"year","lev","lat","lon"/))
filevardef(ncdf,"b1y",typeof(b1y),(/"year","lev","lat","lon"/))
filevardef(ncdf,"b3x",typeof(b3x),(/"year","lev","lat","lon"/))
filevardef(ncdf,"b3y",typeof(b3y),(/"year","lev","lat","lon"/))
filevardef(ncdf,"t3x",typeof(t3x),(/"year","lev","lat","lon"/))
filevardef(ncdf,"t3y",typeof(t3y),(/"year","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"b2",b2)
filevarattdef(ncdf,"b1x",b1x)
filevarattdef(ncdf,"b1y",b1y)
filevarattdef(ncdf,"b3x",b3x)
filevarattdef(ncdf,"b3y",b3y)
filevarattdef(ncdf,"t3x",t3x)
filevarattdef(ncdf,"t3y",t3y)

ncdf->lev  = (/toint(lev)/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->b2 = (/b2/)
ncdf->b1x = (/b1x/)
ncdf->b1y = (/b1y/)
ncdf->b3x = (/b3x/)
ncdf->b3y = (/b3y/)
ncdf->t3x = (/t3x/)
ncdf->t3y = (/t3y/)
end

