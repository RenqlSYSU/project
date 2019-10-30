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
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_daily_eddy_flux.nc"
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
lats = -20 ;-90 ;
latn = 60  ; 90 ;
lonl = 20  ; 0  ;
lonr = 220 ;360 ;
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
month_s = 5
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
t     = new((/nyear,ntime,nlev,nlat,nlon/),float)
w     = new((/nyear,ntime,nlev,nlat,nlon/),float)
Qeddy6 = new((/nyear,4,ntime,nlev,nlat,nlon/),double)
A6     = new((/nyear,2,ntime,nlev,nlat,nlon/),double)
Qeddy  = new((/      4,ntime,nlev,nlat,nlon/),double)
A      = new((/      2,ntime,nlev,nlat,nlon/),double)

iyear = 0
itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
coslat = conform(u,cos(lat),ilat)
sinlat = conform(u,sin(lat),ilat)
lev5d  = conform(u, lev*100, ilev)
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 0,ncase-1,1 
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
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

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
Qeddy6(:,0,:,:,:,:) = - (1.0/a/coslat) * center_finite_diff_n(edd_u*edd_t       ,lon,cyclic,opt,ilon)
Qeddy6(:,1,:,:,:,:) = - (1.0/a/coslat) * center_finite_diff_n(edd_v*edd_t*coslat,lat,cyclic,opt,ilat) 
Qeddy6(:,2,:,:,:,:) = - center_finite_diff_n(edd_w*edd_t,lev*100,cyclic,opt,ilev) 
Qeddy6(:,3,:,:,:,:) =   (R/cp)*(edd_w*edd_t)/lev5d 
Qeddy = dim_avg_n(Qeddy6,iyear)
print(Qeddy(0,15,5,30,30))

vort = center_finite_diff_n(edd_v,lon,cyclic,opt,ilon)- coslat*center_finite_diff_n(edd_u,lat,cyclic,opt,ilat) + sinlat*edd_u 
A6(:,0,:,:,:,:) = -(1.0/a/a/coslat/coslat)*center_finite_diff_n(edd_u*vort, lon, cyclic, opt, ilon)
A6(:,1,:,:,:,:) = -(1.0/a/a/coslat       )*center_finite_diff_n(edd_v*vort, lat, cyclic, opt, ilat)
A = dim_avg_n(A6,iyear)
print(A(0,15,5,30,30))

Qeddy@long_name = "four terms of transient eddy heating, -dutdx, -dvtdy, -dwtdz, wt"
Qeddy@units     = "K/s" 
A@long_name = "two terms of diabatic heating, -duvortdx, -dvvortdy"
A@units     = "s^-2" 

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
dimNames = (/"term1","term2","time","lev","lat","lon"/);
dimSizes = (/ 4     , 2     ,ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False  ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(u&time),getvardims(u&time))
filevardef(ncdf,"lev" ,typeof(u&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
filevardef(ncdf,"Qeddy",typeof(Qeddy),(/"term1","time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A),(/"term2","time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lev",u&level)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"Qeddy",Qeddy)
filevarattdef(ncdf,"A"    ,A)

ncdf->time = (/u&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->Qeddy = (/Qeddy/)
ncdf->A     = (/A/)
end do
end

