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

if(case_option.eq.0) then
    case    = (/"CTRL","NG24h","NG6h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_MJJA_daily_adevection.nc"
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_MJJA_daily_advection.nc"
end if

var_name= (/"U","V","T"/)
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
month_s = 5
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
vars:= f->$var_name(0)$(select_time,{lev},{lats:latn},{lonl:lonr})
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)

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

do nf = 2,2,1 
if((case_option.eq.1).and.(nf.eq.1))
    years  = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
t     = new((/nyear,ntime,nlev,nlat,nlon/),float)

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
end do
clm_u = dim_avg_n(u,0)
clm_v = dim_avg_n(v,0)
clm_t = dim_avg_n(t,0)
delete([/u,v,t/])
printVarSummary(clm_u)
clm_u = linmsg_n(clm_u,-1,ilon)
clm_v = linmsg_n(clm_v,-1,ilon)
clm_t = linmsg_n(clm_t,-1,ilon)

;==========================================================
;calc advection forcing 
;============================================================
coslat = conform(clm_u,cos(lat),ilat)
lev4d  = conform(clm_u, lev*100,ilev)
sigm0 = (R/cp)*(clm_t/lev4d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev)
delete(lev4d)

vort  = (1/a/coslat)*(center_finite_diff_n(clm_v,lon,loncyc,opt,ilon) - center_finite_diff_n(clm_u*coslat,lat,cyclic,opt,ilat))
b1 = - (clm_u/a/coslat)*center_finite_diff_n(vort,lon,loncyc,opt,ilon) \
     - (clm_v/a       )*center_finite_diff_n(vort,lat,cyclic,opt,ilat)
print(vort(1,10,20:30))
print(b1(0,10,20:30))
delete(vort)

b2 = - 2*(2.0*pi/24.0/3600.0)*coslat*clm_v/a

f0   = conform(clm_u,2*(2*pi/24.0/3600.0)*sin(lat),ilat)
advT = (clm_u/a/coslat)*center_finite_diff_n(clm_t,lon,loncyc,opt,ilon) + \
       (clm_v/a       )*center_finite_diff_n(clm_t,lat,cyclic,opt,ilat)
delete([/clm_u,clm_v,clm_t,coslat/])
b3 = f0 * center_finite_diff_n(advT/sigm0,lev*100,cyclic,opt,ilev)
delete([/advT,f0/])

sigm0@long_name = "one kind of static stability measured by temperature"
sigm0@units     = "K/Pa" 
b1@long_name = "Horizontal advection of relative vorticity"
b1@units     = "s^-2" 
b2@long_name = "Horizontal advection of absolute vorticity"
b2@units     = "s^-2" 
b3@long_name = "Vertical gradient of temperature horizontal advection"
b3@units     = "s^-2" 

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
dimSizes = (/ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"sigm0",typeof(sigm0),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b1",typeof(b1),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b2",typeof(b2),(/"time","lev","lat","lon"/))
filevardef(ncdf,"b3",typeof(b3),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm0",sigm0)
filevarattdef(ncdf,"b1",b1)
filevarattdef(ncdf,"b2",b2)
filevarattdef(ncdf,"b3",b3)

ncdf->time = (/vars&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm0 = (/sigm0/)
ncdf->b1 = (/b1/)
ncdf->b2 = (/b2/)
ncdf->b3 = (/b3/)
end do
end

