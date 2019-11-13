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
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_eddy_flux2.nc"
    else
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_eddy_flux_" + cb + "-" + ca + "day.nc"
    end if  
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux.nc"
    else                                                                                              
        fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_eddy_flux_" + cb + "-" + ca + "day.nc"
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
month_s = 1
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
nmonth= 4
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
Qeddy_m  = new((/   3,nmonth,nlev,nlat,nlon/),float)
A_m      = new((/   3,nmonth,nlev,nlat,nlon/),float)

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

do nf = 0,ncase-1,2
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
delete([/clm_u,clm_v,clm_t,clm_w/])

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
Qeddy6 = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
Qeddy6(:,0,:,:,:,:) = edd_u*edd_t
Qeddy6(:,1,:,:,:,:) = edd_v*edd_t
Qeddy6(:,2,:,:,:,:) = edd_w*edd_t
Qeddy = tofloat(dim_avg_n(Qeddy6,iyear))
delete(Qeddy6)
print(Qeddy(0,15,5,30,30))

A6   = new((/nyear,3,ntime,nlev,nlat,nlon/),double)
vort = center_finite_diff_n(edd_v,lon,cyclic,opt,ilon) - coslat*center_finite_diff_n(edd_u,lat,cyclic,opt,ilat) + sinlat*edd_u 
A6(:,0,:,:,:,:) = edd_u*vort
A6(:,1,:,:,:,:) = edd_v*vort
A6(:,2,:,:,:,:) = edd_w*vort
A = tofloat(dim_avg_n(A6,iyear))
delete([/A6,vort/])
delete([/edd_u,edd_v,edd_t,edd_w/])
print(A(0,15,5,30,30))

Qeddy@long_name = "three terms of transient eddy heating, ut, vt, wt"
Qeddy@units     = "K*m/s" 
A@long_name = "three terms of diabatic heating, uvort, vvort, wvort"
A@units     = "m*s^-2" 

time1 = cd_calendar(f->time(select_time),0)
do nm = 0, 3, 1
    month1_s = 2*nm + 1
    month1_e = 2*nm + 2
    select_time1 := ind(time1(:,1).ge.month1_s.and.time1(:,1).le.month1_e)
    Qeddy_m(:,nm,:,:,:) = dim_avg_n_Wrap(Qeddy(:,select_time1,:,:,:),1)
        A_m(:,nm,:,:,:) = dim_avg_n_Wrap(    A(:,select_time1,:,:,:),1)
end do 
delete([/Qeddy,A/])

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
dimNames = (/"term1","time","lev","lat","lon"/);
dimSizes = (/ 3     ,nmonth,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"Qeddy",typeof(Qeddy_m),(/"term1","time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A_m)    ,(/"term1","time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"Qeddy",Qeddy_m)
filevarattdef(ncdf,"A"    ,A_m)

ncdf->time = (/12,34,56,78/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->Qeddy = (/Qeddy_m/)
ncdf->A     = (/A_m/)
end do
end

