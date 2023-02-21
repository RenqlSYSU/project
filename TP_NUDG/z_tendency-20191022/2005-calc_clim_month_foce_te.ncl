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
nmonth = 2
season = False;True ;
geo = True ;False ;if true, only can calc 10N 

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_forc_tei_season.nc"
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_forc_te"
end if

if(case_option.eq.2) then
    case    = (/"F2000","NGTR","NUDG","NG58"/)
    filein  = path +(/"/model_data/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model_data/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    if (season) then 
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_season_forc_te.nc"
    else
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_2month_forc_te.nc"
    end if
end if

var_name = (/"U","V","T","OMEGA","Z3"/)
ncase    = dimsizes(case)
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
zlev = (/1000,925,850/) ;hPa ,700,600,500,400
lev2 = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;hPa 13levels
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
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
nzlev = dimsizes(zlev)
lev := lev3
nlev  = dimsizes(lev)
Qeddl_m    = new((/nmonth,nlev,nlat,nlon/),double)
Qeddh_m    = new((/nmonth,nlev,nlat,nlon/),double)
Al_m       = new((/nmonth,nlev,nlat,nlon/),double)
Ah_m       = new((/nmonth,nlev,nlat,nlon/),double)
sigm0_m    = new((/nmonth,nlev,nlat,nlon/),double)

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

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years  = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nm = 0, nmonth-1, 1
if(season) then 
if(nm.eq.0) then 
month_s = 12  
month_e = 2
select_time := ind(time(:,1).eq.12.or.time(:,1).le.2)
else
month_s = 3*nm ;+ 1
month_e = 3*nm + 2
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
else
month_s = 2*nm + 5
month_e = 2*nm + 6
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
w     := new((/nyear,ntime,nlev ,nlat,nlon/),float)
z     := new((/nyear,ntime,nzlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
print("handle with " + case(nf) + " " + years(nt))
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
w(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 4
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f      = addfile(files(nt),"r")
z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{zlev},{lats:latn},{lonl:lonr})
end do
u = linmsg_n(u,-1,ilon)
v = linmsg_n(v,-1,ilon)
t = linmsg_n(t,-1,ilon)
w = linmsg_n(w,-1,ilon)
z = linmsg_n(z,-1,ilon)
printVarSummary(w)

if(geo) then
coslat = conform(z(:,:,0,:,:),cos(lat),ilat-1)
f0     = conform(z(:,:,0,:,:),2*(2*pi/24.0/3600.0)*sin(lat),ilat-1)
z = z*g
lat1 = vars&lat
lat_ind = ind(lat1.le.90.0)
lat2 = lat(lat_ind)
printVarSummary(lat_ind)
do nlv = 0,1,1
    nl = 1 - nlv
    z1 = z(:,:,nl+1,lat_ind,:) + R*t(:,:,nl,lat_ind,:)*(lev(nl+1)-lev(nl))*100/lev(nl)/100
    u1 = -center_finite_diff_n(z1,lat2,cyclic,opt,ilat-1)/f0(:,:,lat_ind,:)/a        ;(:,:,nl,lat_ind,:)
    v1 =  center_finite_diff_n(z1,lon ,loncyc,opt,ilon-1)/f0(:,:,lat_ind,:)/a/coslat(:,:,lat_ind,:) ;(:,:,nl,lat_ind,:)
    t1 = -lev(nl)*100*(z(:,:,nl+1,lat_ind,:)-z1)/((lev(nl+1)-lev(nl))*100)/R
    print(max(abs(z(:,:,nl,lat_ind,:)-z1)))
    print(max(abs(u(:,:,nl,lat_ind,:)-u1)))
    print(max(abs(v(:,:,nl,lat_ind,:)-v1)))
    print(max(abs(t(:,:,nl,lat_ind,:)-t1)))
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

lev := lev3
nlev = dimsizes(lev)
Qeddy := new((/nyear,ntime,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(u    ,lev2,ilev)
p_out  = conform(Qeddy,lev3,ilev)
u := int2p_n(p_in,u,p_out,linlog,ilev)
v := int2p_n(p_in,v,p_out,linlog,ilev)
t := int2p_n(p_in,t,p_out,linlog,ilev)
w := int2p_n(p_in,w,p_out,linlog,ilev)
delete([/p_in,p_out,linlog/])

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
sigm0_m(nm,:,:,:) = dim_avg_n(sigm0,itime-1)
delete([/lev4d,clm_t,sigm0/])

coslat = conform(hgh_u,cos(lat),ilat)
tanlat = conform(hgh_u,tan(lat),ilat)
lev5d  = conform(hgh_u, lev*100,ilev)

Qeddh := - (1.0/a/coslat) * center_finite_diff_n(hgh_u*hgh_t       ,lon,loncyc,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(hgh_v*hgh_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(hgh_w*hgh_t)/lev5d - center_finite_diff_n(hgh_w*hgh_t,lev*100,cyclic,opt,ilev) 
Qeddh := dim_avg_n(Qeddh,itime)
Qeddh_m(nm,:,:,:) = dim_avg_n(Qeddh,iyear)
delete(Qeddh)

Qeddl := - (1.0/a/coslat) * center_finite_diff_n(low_u*low_t       ,lon,loncyc,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(low_v*low_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(low_w*low_t)/lev5d - center_finite_diff_n(low_w*low_t,lev*100,cyclic,opt,ilev) 
Qeddl := dim_avg_n(Qeddl,itime)
Qeddl_m(nm,:,:,:) = dim_avg_n(Qeddl,iyear)
printVarSummary(Qeddl)
print(Qeddl(0,5,30,30))
delete([/lev5d,Qeddl/])

Ah:= -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_v              ,lon    ,loncyc,opt,ilon),lon    ,loncyc,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_v*hgh_v*coslat       ,lon    ,loncyc,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(hgh_v*hgh_w              ,lon    ,loncyc,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(hgh_u*hgh_u              ,lon    ,loncyc,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_u              ,lat    ,cyclic,opt,ilat),lon    ,loncyc,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(hgh_u*hgh_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(hgh_u*hgh_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     -(tanlat/a                )*                     center_finite_diff_n(hgh_u*hgh_w              ,lev*100,cyclic,opt,ilev)
Ah := dim_avg_n(Ah,itime)
Ah_m(nm,:,:,:) = dim_avg_n(Ah,iyear)
printVarSummary(Ah)
print(Ah(0,5,30,30))
delete([/hgh_u,hgh_v,hgh_t,hgh_w,Ah/])

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
Al_m(nm,:,:,:) = dim_avg_n(Al,iyear)
printVarSummary(Al)
print(Al(0,5,30,30))
delete([/low_u,low_v,low_t,low_w,Al/])
delete([/coslat,tanlat/])
end do

Qeddh_m@long_name = "Climatology 2-10-day transient eddy heating"
Qeddh_m@units     = "K/s" 
Qeddl_m@long_name = "Climatology 10-90-day transient eddy heating"
Qeddl_m@units     = "K/s" 
Ah_m@long_name = "Climatology 2-10-day diabatic heating"
Ah_m@units     = "s^-2" 
Al_m@long_name = "Climatology 10-90-day diabatic heating"
Al_m@units     = "s^-2" 

;==========================================================
;calc Qd foring and Qeddy forcing, dimension is 4d 
;============================================================
lev4d = conform(Ah_m, lev*100,ilev-1)
sigm1 = sigm0_m*R/lev4d
sigm1_m = conform(lev4d, wgt_areaave(sigm1,1.0,1.0,0), (/0,ilev-1/))
f0    = conform(Ah_m, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-1)
f_Qeddh  = -f0*R*center_finite_diff_n((Qeddh_m/sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
f_Qeddl  = -f0*R*center_finite_diff_n((Qeddl_m/sigm1_m/lev4d),lev*100,cyclic,opt,ilev-1)
delete(lev4d)
delete(f0)
delete(sigm1_m)

f_Qeddh@long_name = "Climatology 2-10-day transient eddy heating"
f_Qeddh@units     = "s^-2"                 
f_Qeddl@long_name = "Climatology 10-90-day transient eddy heating"
f_Qeddl@units     = "s^-2" 
sigm1@long_name = "One kind of static stability measured by potential temperature"
sigm1@units     = "J/(Pa^2*Kg)" 
print(lat(0:9))

;========================================================================

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
dimSizes = (/nmonth,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"sigm1",typeof(sigm1),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Al"   ,typeof(Al_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Ah"   ,typeof(Ah_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qeddh",typeof(Qeddh_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qeddl",typeof(Qeddl_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddh",typeof(f_Qeddh),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddl",typeof(f_Qeddl),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm1",sigm1)
filevarattdef(ncdf,"Al"   ,Al_m)
filevarattdef(ncdf,"Ah"   ,Ah_m)
filevarattdef(ncdf,"Qeddl",Qeddl_m)
filevarattdef(ncdf,"Qeddh",Qeddh_m)
filevarattdef(ncdf,"f_Qeddl",f_Qeddl)
filevarattdef(ncdf,"f_Qeddh",f_Qeddh)

ncdf->lev  = (/toint(lev)/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm1 = (/sigm1/)
ncdf->Al    = (/Al_m/)
ncdf->Ah    = (/Ah_m/)
ncdf->Qeddl = (/Qeddl_m/)
ncdf->Qeddh = (/Qeddh_m/)
ncdf->f_Qeddl = (/f_Qeddl/)
ncdf->f_Qeddh = (/f_Qeddh/)
end do
end

