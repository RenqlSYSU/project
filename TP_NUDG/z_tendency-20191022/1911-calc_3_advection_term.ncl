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
nmonth = 4
season = False 
geo = True ;False ;if true, only can calc 10N 

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_advect2.nc"
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_advect.nc"
end if

if(case_option.eq.2) then
    case    = (/"F2000","NUDG","NG58"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_advect.nc"
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
zlev = (/1000,925,850/) ;hPa ,700,600,500,400
lev2 = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;13 levels,hPa
lev3 = fspan(1000,100,19)
lev = lev2

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
nzlev = dimsizes(zlev)
lev := lev3
nlev  = dimsizes(lev)
b1_m  = new((/nmonth,nlev,nlat,nlon/),float)
b2_m  = new((/nmonth,nlev,nlat,nlon/),float)
b3_m  = new((/nmonth,nlev,nlat,nlon/),float)
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
month_s = 2*nm + 1
month_e = 2*nm + 2
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
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

;nv = 3
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
;z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
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
b1 := new((/ntime,nlev,nlat,nlon/),double)

linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
p_in   = conform(clm_u,lev2,ilev)
p_out  = conform(b1,lev3,ilev)
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
b1 = - (clm_u/a/coslat)*center_finite_diff_n(vort,lon,loncyc,opt,ilon) 
b2 = - (clm_v/a       )*center_finite_diff_n(vort,lat,cyclic,opt,ilat)-2*(2.0*pi/24.0/3600.0)*coslat*clm_v/a
print(vort(1,9,10,20:30))
print(b1(1,9,10,20:30))
vort_m(nm,:,:,:)  = tofloat(dim_avg_n(vort,0))
b1_m(nm,:,:,:)    = tofloat(dim_avg_n(b1,0))
b2_m(nm,:,:,:)    = tofloat(dim_avg_n(b2,0))
sigm0_m(nm,:,:,:) = tofloat(dim_avg_n(sigm0,0))
delete([/vort,b1,b2/])

sigm0 = conform(clm_u,wgt_areaave(sigm0,1.0,1.0,0),(/itime,ilev/))
f0   = conform(clm_u,2*(2*pi/24.0/3600.0)*sin(lat),ilat)
advT = (clm_u/a/coslat)*center_finite_diff_n(clm_t,lon,loncyc,opt,ilon) + \
       (clm_v/a       )*center_finite_diff_n(clm_t,lat,cyclic,opt,ilat)
delete([/clm_u,clm_v,clm_t,coslat/])
b3 = f0 * center_finite_diff_n(advT/sigm0,lev*100,cyclic,opt,ilev)
b3_m(nm,:,:,:) = tofloat(dim_avg_n(b3,0))
delete([/advT,f0,b3,sigm0/])
end do

sigm0_m@long_name = "one kind of static stability measured by temperature"
sigm0_m@units     = "K/Pa" 
vort_m@long_name = "one kind of static stability measured by temperature"
vort_m@units     = "K/Pa" 
b1_m@long_name = "x transport of Absolute vorticity"
b1_m@units     = "s^-2" 
b2_m@long_name = "y transport of Absolute vorticity"
b2_m@units     = "s^-2" 
b3_m@long_name = "Vertical gradient of temperature horizontal advection"
b3_m@units     = "s^-2" 

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
filevardef(ncdf,"b1",typeof(b1_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b2",typeof(b2_m),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b3",typeof(b3_m),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm0",sigm0_m)
filevarattdef(ncdf,"vort" ,vort_m)
filevarattdef(ncdf,"b1",b1_m)
filevarattdef(ncdf,"b2",b2_m)
filevarattdef(ncdf,"b3",b3_m)

ncdf->lev  = (/toint(lev)/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm0 = (/sigm0_m/)
ncdf->vort  = (/vort_m/)
ncdf->b1 = (/b1_m/)
ncdf->b2 = (/b2_m/)
ncdf->b3 = (/b3_m/)
end do
end

