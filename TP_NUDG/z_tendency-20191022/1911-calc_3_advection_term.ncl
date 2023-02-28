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
path    = "/home/ys17-23/Extension2/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment ; 2 means F2000
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
nmonth = 3
season = True ;False 
geo = False;True ;if true, only can calc 10N 

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_advect2_season.nc"
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_advect.nc"
end if

if(case_option.eq.2) then
    case    = (/"F2000","NGTR","NG58"/)
    filein  = path +(/"/model_data/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model_data/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_season_advect2.nc"
end if

var_name= (/"U","V","T","Z3"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats =  10 ;-90 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;
lev2 = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;13 levels,hPa
lev3 = fspan(1000,100,19)
lev = lev2

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
time = cd_calendar(f->time,0)
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
nzlev = 5 
lev := lev3
nlev  = dimsizes(lev)
b1x_m  = new((/nmonth,nlev,nlat,nlon/),float)
b1y_m  = new((/nmonth,nlev,nlat,nlon/),float)
b2_m   = new((/nmonth,nlev,nlat,nlon/),float)
t3x_m  = new((/nmonth,nlev,nlat,nlon/),float)
t3y_m  = new((/nmonth,nlev,nlat,nlon/),float)
b3x_m = new((/nmonth,nlev,nlat,nlon/),float)
b3y_m = new((/nmonth,nlev,nlat,nlon/),float)
sigm0_m  = new((/nmonth,nlev,nlat,nlon/),float)
vort_m   = new((/nmonth,nlev,nlat,nlon/),float)

itime = 0
ilev  = 1
ilat  = 2
ilon  = 3
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
loncyc = True ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years  = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nm = 0,nmonth-1, 1
if(nm.eq.0.and.season) then 
month_s = 12  
month_e = 2
select_time := ind(time(:,1).eq.12.or.time(:,1).le.2)
else
month_s = 3*nm
month_e = 3*nm + 2
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
ntime = dimsizes(select_time)
lev := lev2
nlev  = dimsizes(lev)
print("===================================================================")
print("month from " + month_s + " to " + month_e + ", ntime = " + ntime)
u     := new((/nyear,ntime,nlev ,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev ,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev ,nlat,nlon/),float)
z     := new((/nyear,ntime,nzlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("handle with " + case(nf) + " " + years(nt))
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev(0:(nzlev-1))},{lats:latn},{lonl:lonr})
end do

if(geo) then
coslat = conform(z(:,:,0,:,:),cos(lat),ilat)
f0     = conform(z(:,:,0,:,:),2*(2*pi/24.0/3600.0)*sin(lat),ilat)
z = z*g
lat1 = vars&lat
lat_ind = ind(lat1.le.90.0)
lat2 = lat(lat_ind)
do nlv = 0,1,1
    nl = 1 - nlv
    z1 = z(:,:,nl+1,lat_ind,:) + R*t(:,:,nl,lat_ind,:)*(lev(nl+1)-lev(nl))*100/lev(nl)/100
    u1 = -center_finite_diff_n(z1,lat2,cyclic,opt,ilat)/f0(:,:,lat_ind,:)/a        ;(:,:,nl,lat_ind,:)
    v1 =  center_finite_diff_n(z1,lon ,loncyc,opt,ilon)/f0(:,:,lat_ind,:)/a/coslat(:,:,lat_ind,:) ;(:,:,nl,lat_ind,:)
    t1 = -lev(nl)*100*(z(:,:,nl+1,lat_ind,:)-z1)/((lev(nl+1)-lev(nl))*100)/R
    print(max(abs(z(:,:,nl,:,:)-z1)))
    print(max(abs(u(:,:,nl,:,:)-u1)))
    print(max(abs(v(:,:,nl,:,:)-v1)))
    print(max(abs(t(:,:,nl,:,:)-t1)))
    print(z(1,1,nl,10,10:20)-z1(1,1,10,10:20))
    print(u(1,1,nl,10,10:20)-u1(1,1,10,10:20))
    print(v(1,1,nl,10,10:20)-v1(1,1,10,10:20))
    print(t(1,1,nl,10,10:20)-t1(1,1,10,10:20))
    z(:,:,nl,lat_ind,:) = (/tofloat(z1)/)
    u(:,:,nl,lat_ind,:) = (/tofloat(u1)/)
    v(:,:,nl,lat_ind,:) = (/tofloat(v1)/)
    ;t(:,:,nl,lat_ind,:) = (/tofloat(t1)/)
end do
delete([/coslat,f0,u1,v1,z1,t1,z,lat1,lat_ind,lat2/])
end if
printVarSummary(u)

clm_u = dim_avg_n(u,0)
clm_v = dim_avg_n(v,0)
clm_t = dim_avg_n(t,0)
delete([/u,v,t/])
printVarSummary(clm_u)
clm_u = linmsg_n(clm_u,-1,ilon)
clm_v = linmsg_n(clm_v,-1,ilon)
clm_t = linmsg_n(clm_t,-1,ilon)
print(clm_u(1,9,10,20:30))
print(clm_v(1,9,10,20:30))
print(clm_t(1,9,10,20:30))

lev := lev3
nlev = dimsizes(lev)
b1x := new((/ntime,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(clm_u,lev2,ilev)
p_out  = conform(b1x,lev3,ilev)
clm_u := int2p_n(p_in,clm_u,p_out,linlog,ilev)
clm_v := int2p_n(p_in,clm_v,p_out,linlog,ilev)
clm_t := int2p_n(p_in,clm_t,p_out,linlog,ilev)
delete([/p_in,p_out,linlog/])

;==========================================================
;calc advection forcing 
;============================================================
coslat = conform(clm_u,cos(lat),ilat)
lev4d  = conform(clm_u, lev*100,ilev)
sigm0 = (R/cp)*(clm_t/lev4d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev)
delete(lev4d)

;vort = (1/a/coslat)*(center_finite_diff_n(clm_v,lon,loncyc,opt,ilon) - center_finite_diff_n(clm_u*coslat,lat,cyclic,opt,ilat))
vort = uv2vr_cfd(clm_u,clm_v,vars&lat,vars&lon,3) 
b1x = - (clm_u/a/coslat)*center_finite_diff_n(vort,lon,loncyc,opt,ilon) 
b1y = - (clm_v/a       )*center_finite_diff_n(vort,lat,cyclic,opt,ilat)
b2 = -2*(2.0*pi/24.0/3600.0)*coslat*clm_v/a
print(vort(1,9,10,20:30))
print(b1x(1,9,10,20:30))
vort_m(nm,:,:,:)  = tofloat(dim_avg_n(vort,0))
b1x_m(nm,:,:,:)   = tofloat(dim_avg_n(b1x,0))
b1y_m(nm,:,:,:)   = tofloat(dim_avg_n(b1y,0))
b2_m(nm,:,:,:)    = tofloat(dim_avg_n(b2,0))
sigm0_m(nm,:,:,:) = tofloat(dim_avg_n(sigm0,0))
delete([/vort,b1x,b1y,b2/])

sigm0 = conform(clm_u,wgt_areaave(sigm0,1.0,1.0,0),(/itime,ilev/))
f0   = conform(clm_u,2*(2*pi/24.0/3600.0)*sin(lat),ilat)
t3x  = -1*(clm_u/a/coslat)*center_finite_diff_n(clm_t,lon,loncyc,opt,ilon)
t3y  = -1*(clm_v/a       )*center_finite_diff_n(clm_t,lat,cyclic,opt,ilat)
delete([/clm_u,clm_v,clm_t,coslat/])
b3x = -1*f0*center_finite_diff_n(t3x/sigm0,lev*100,cyclic,opt,ilev)
b3y = -1*f0*center_finite_diff_n(t3y/sigm0,lev*100,cyclic,opt,ilev)
b3x_m(nm,:,:,:) = tofloat(dim_avg_n(b3x,0))
b3y_m(nm,:,:,:) = tofloat(dim_avg_n(b3y,0))
t3x_m(nm,:,:,:) = tofloat(dim_avg_n(t3x,0))
t3y_m(nm,:,:,:) = tofloat(dim_avg_n(t3y,0))
delete([/f0,b3x,b3y,t3x,t3y,sigm0/])
end do

sigm0_m@long_name = "one kind of static stability measured by temperature"
sigm0_m@units     = "K/Pa" 
vort_m@long_name = "one kind of static stability measured by temperature"
vort_m@units     = "K/Pa" 
b1x_m@long_name = "x transport of Absolute vorticity"
b1x_m@units     = "s^-2" 
b2_m@long_name = "y transport of Absolute vorticity"
b2_m@units     = "s^-2" 
b3x_m@long_name = "Vertical gradient of temperature horizontal advection"
b3x_m@units     = "s^-2" 

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
fileatt@discription   = "months from " + month_s + " to " + month_e + "daily three advection terms"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lev","lat","lon"/);
dimSizes = (/nmonth,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"sigm0",typeof(sigm0_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"vort" ,typeof(vort_m) ,(/"time","lev","lat","lon"/))
filevardef(ncdf,"b1x",typeof(b1x_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b1y",typeof(b1y_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b2",typeof(b2_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b3x",typeof(b3x_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b3y",typeof(b3y_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"t3x",typeof(t3x_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"t3y",typeof(t3y_m),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm0",sigm0_m)
filevarattdef(ncdf,"vort" ,vort_m)
filevarattdef(ncdf,"b1x",b1x_m)
filevarattdef(ncdf,"b1y",b1y_m)
filevarattdef(ncdf,"b2",b2_m)
filevarattdef(ncdf,"b3x",b3x_m)
filevarattdef(ncdf,"b3y",b3y_m)
filevarattdef(ncdf,"t3x",t3x_m)
filevarattdef(ncdf,"t3y",t3y_m)

ncdf->lev  = (/toint(lev)/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm0 = (/sigm0_m/)
ncdf->vort  = (/vort_m/)
ncdf->b1x = (/b1x_m/)
ncdf->b1y = (/b1y_m/)
ncdf->b2 = (/b2_m/)
ncdf->b3x = (/b3x_m/)
ncdf->b3y = (/b3y_m/)
ncdf->t3x = (/t3x_m/)
ncdf->t3y = (/t3y_m/)
end do
end

