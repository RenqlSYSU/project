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
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_daily_Qd_Qeddy_A-global.nc"
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
lats = -20;-90 ;
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
Qd    = new((/nyear,ntime,nlev,nlat,nlon/),float)
Qeddy = new((/nyear,ntime,nlev,nlat,nlon/),float)
A     = new((/nyear,ntime,nlev,nlat,nlon/),float)
f_Qd    = new((/    ntime,nlev,nlat,nlon/),double)
f_Qeddy = new((/    ntime,nlev,nlat,nlon/),double)

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

do nf = 0,0,1 
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
t     = new((/nyear,ntime,nlev,nlat,nlon/),float)
w     = new((/nyear,ntime,nlev,nlat,nlon/),float)

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

coslat = conform(Qd,cos(lat),ilat)
tanlat = conform(Qd,tan(lat),ilat)
lev5d  = conform(Qd, lev*100, ilev)

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
sigm0 = (R/cp)*(clm_t/lev5d) - center_finite_diff_n(clm_t,lev*100,cyclic,opt,ilev)
Qd   :=   (clm_u/a/coslat) * center_finite_diff_n(clm_t,lon,cyclic,opt,ilon) \
        + (clm_v/a       ) * center_finite_diff_n(clm_t,lat,cyclic,opt,ilat) \
        - sigm0 * clm_w 
delete([/clm_u,clm_v,clm_t,clm_w/])

Qeddy := - (1.0/a/coslat) * center_finite_diff_n(edd_u*edd_t       ,lon,cyclic,opt,ilon) \
         - (1.0/a/coslat) * center_finite_diff_n(edd_v*edd_t*coslat,lat,cyclic,opt,ilat) \
         + (R/cp)*(edd_w*edd_t)/lev5d - center_finite_diff_n(edd_w*edd_t,lev*100,cyclic,opt,ilev) 
Qeddy := dim_avg_n(Qeddy,iyear)
print(Qeddy(15,5,30,30))

Qd   := dim_avg_n(Qd   ,iyear) - Qeddy
print(Qd  (15,5,30,30))

A := -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v              ,lon    ,cyclic,opt,ilon),lon    ,cyclic,opt,ilon)\
     -(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_v*edd_v*coslat       ,lon    ,cyclic,opt,ilon),lat    ,cyclic,opt,ilat)\
     -(1.0   /a  /coslat       )*center_finite_diff_n(center_finite_diff_n(edd_v*edd_w              ,lon    ,cyclic,opt,ilon),lev*100,cyclic,opt,ilev)\
     -(tanlat/a/a/coslat       )*                     center_finite_diff_n(edd_u*edd_u              ,lon    ,cyclic,opt,ilon)                         \
     +(1.0   /a/a/coslat       )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_u              ,lat    ,cyclic,opt,ilat),lon    ,cyclic,opt,ilon)\
     +(1.0   /a/a/coslat/coslat)*center_finite_diff_n(center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat),lat    ,cyclic,opt,ilat)\
     +(tanlat/a/a/coslat/coslat)*                     center_finite_diff_n(edd_u*edd_v*coslat*coslat,lat    ,cyclic,opt,ilat)                         \
     +(1.0   /a                )*center_finite_diff_n(center_finite_diff_n(edd_u*edd_w              ,lat    ,cyclic,opt,ilat),lev*100,cyclic,opt,ilev)\
     -(tanlat/a                )*                     center_finite_diff_n(edd_u*edd_w              ,lev*100,cyclic,opt,ilev)
A := dim_avg_n(A,iyear)
print(A(15,5,30,30))
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
sigm1 = dim_avg_n(sigm0*R/lev5d,iyear)
delete(sigm0)
lev4d = dim_avg_n(lev5d,iyear)
delete(lev5d)
f0    = conform(sigm1, 2*(2*pi/24.0/3600.0)*sin(lat), ilat-1)

;do ilat = 0, nlat-1
;if (abs(lat(ilat)*180.0/pi ).lt. 10. ) then
;f0(:,:,ilat,:)= A@_FillValue
;end if
;end do

f_Qd     = -f0*R*center_finite_diff_n((Qd   /sigm1/lev4d),lev*100,cyclic,opt,ilev-1)
f_Qeddy  = -f0*R*center_finite_diff_n((Qeddy/sigm1/lev4d),lev*100,cyclic,opt,ilev-1)
delete(lev4d)

f_Qeddy@long_name = "Climatology daily transient eddy heating"
f_Qeddy@units     = "s^-2" 
f_Qd@long_name    = "Climatology daily diabatic heatinig"
f_Qd@units        = "s^-2" 
sigm1@long_name = "One kind of static stability measured by potential temperature"
sigm1@units     = "J/(Pa^2*Kg)" 
print(f_Qeddy(15,5,30,30))
print(f_Qd   (15,5,30,30))
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
fileatt@discription   = "months from " + month_s + " to " + month_e + "daily three forcing factors, Qd, Qeddy, A"
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
filevardef(ncdf,"sigm1",typeof(sigm1),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qeddy",typeof(Qeddy),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qd"   ,typeof(Qd),(/"time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qd",typeof(f_Qd),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddy",typeof(f_Qeddy),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"sigm1",sigm1)
filevarattdef(ncdf,"Qeddy",Qeddy)
filevarattdef(ncdf,"Qd"   ,Qd)
filevarattdef(ncdf,"A"    ,A)
filevarattdef(ncdf,"f_Qeddy",f_Qeddy)
filevarattdef(ncdf,"f_Qd"   ,f_Qd)

ncdf->time = (/vars&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->sigm1 = (/sigm1/)
ncdf->Qeddy = (/Qeddy/)
ncdf->Qd    = (/Qd/)
ncdf->A     = (/A/)
ncdf->f_Qeddy = (/f_Qeddy/)
ncdf->f_Qd    = (/f_Qd/)
end do
end

