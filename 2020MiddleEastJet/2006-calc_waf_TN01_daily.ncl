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
var_name = (/"uwnd","vwnd","air","hgt"/)
fileout  = path2 + "ERA-Interim_DJF_TN01_2-year.nc"
year_option = True;False 

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
sclhgt = 7991 ;m 
pi = atan(1.0)*4
print(pi)
years    = ispan(1979,1981,1)
years@calendar = "standard"

;area: east Asia
lats = -90  ;-20
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
nlev  = dimsizes(lev)
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
z     := new((/nyear,ntime,nlev,nlat,nlon/),float)

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
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end if
end do

;the Brunt-Vaisala frequency of atmosphere
theta = pot_temp(lev*100, t, ilev, False)
delete(t)
hgt = z/9.8 ;convert unit from m2/s2 to m
NN  = brunt_vaisala_atm(theta, hgt, 0, ilev)
delete([/theta,hgt/])
NN = NN*abs(NN)
NN = where(abs(NN).lt.0.0000000001, 0.000001, NN)

;noz    = -sclhgt*log(lev/1000.0)

if (any(ismissing(z))) then
    print("hgt has missing")
end if
if (any(ismissing(u))) then
    print("uwnd has missing")
end if
if (any(ismissing(v))) then
    print("vwnd has missing")
end if
z = linmsg_n(z,-1,ilon)

;stream function
f0 = conform(z,2.0*(2.0*pi/24.0/3600.0)*sin(lat),ilat)
sf = u
vp = v
uv2sfvpf ( u, v, sf, vp )
delete(vp)
;sf = z/f0
NN = f0*f0/NN
delete([/z,f0/])

if (any(ismissing(sf))) then
    print("stream function has missing")
end if

;--------------------------------------
;DJF
;---------------------------------------------------
u2      := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2      := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
sf2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
NN2 := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    sf2(ny,0:30,:,:,:) = sf(ny,59:89,:,:,:)
    NN2(ny,0:30,:,:,:) = NN(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    sf2(ny,31:89,:,:,:) = sf(ny+1,0:58,:,:,:)
    NN2(ny,31:89,:,:,:) = NN(ny+1,0:58,:,:,:)
end do
u  := u2
v  := v2
sf := sf2    
NN := NN2
printVarSummary(u)

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
ca = 10.0 ;day
fca = 1.0/ca
fcf = fspan(0, 0.5, ntime/2) 
hgh_cf = ind(fcf.ge.fca)
low_cf = ind(fcf.lt.fca)

wafx  := new((/3,nyear-1,nlev,nlat,nlon/),float)
wafy  := new((/3,nyear-1,nlev,nlat,nlon/),float)
wafz  := new((/3,nyear-1,nlev,nlat,nlon/),float)

do neddy = 0,2,1
if(neddy.eq.0) then
    clm_sf = dim_avg_n(sf,itime)
    edd_sf = sf - conform(sf, clm_sf, (/0,2,3,4/))
    delete([/clm_sf/])
end if
if(neddy.eq.1) then ;low-frequency
    cf = ezfftf_n(sf,itime)
    cf(:,:,hgh_cf,:,:,:) = 0.0  
    edd_sf = ezfftb_n(cf,0,itime+1)
    delete([/cf/])
end if
if(neddy.eq.2) then ;hgh-frequency
    cf = ezfftf_n(sf,itime)
    cf(:,:,low_cf,:,:,:) = 0.0  
    edd_sf = ezfftb_n(cf,0,itime+1)
    delete([/cf/])
end if

dsfdx = center_finite_diff_n(edd_sf,lon,loncyc,opt,ilon)
dsfdy = center_finite_diff_n(edd_sf,lat,cyclic,opt,ilat)

coslat = conform(u,cos(lat),ilat)
coef   = (coslat*conform(u,lev,ilev)/1000.0)/2.0/sqrt(u^2 + v^2)
fx  = coef*((u/a/a/coslat/coslat)*(dsfdx*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lon,loncyc,opt,ilon)) \
           -(v/a/a/coslat       )*(dsfdy*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lat,cyclic,opt,ilat)))
wafx(neddy,:,:,:,:) = dim_avg_n(fx,itime)
delete(fx)

fy  = coef*((u/a/a/coslat)*(dsfdx*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lon,loncyc,opt,ilon)) \
           -(v/a/a       )*(dsfdy*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lat,cyclic,opt,ilat)))
wafy(neddy,:,:,:,:) = dim_avg_n(fy,itime)
delete(fy)

dsfdz = center_finite_diff_n(edd_sf,lev*100,cyclic,opt,ilev)
fz  = coef*NN*((u/a/coslat)*(dsfdx*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lon,loncyc,opt,ilon)) \
              -(v/a       )*(dsfdy*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lat,cyclic,opt,ilat)))
wafz(neddy,:,:,:,:) = dim_avg_n(fz,itime)
delete([/coslat,coef,dsfdx,dsfdy,dsfdz,fz/])
end do

wafx@long_name = "x wave activity flux"
wafx@units   = "m2/s2" 
wafy@long_name = "y wave activity flux"
wafy@units     = "m2/s2" 
wafz@long_name = "z wave activity flux"
wafz@units     = "m*Pa/s2" 

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
fileatt@discription   = "DJF mean wave activity flux based TN01 for total, low-frequency and high-frequency"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"nvar","year" ,"lev","lat","lon"/);
dimSizes = (/  3   ,nyear-1,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(year_option) then 
filevardef(ncdf,"wafx",typeof(wafx),(/"nvar","year","lev","lat","lon"/))
filevardef(ncdf,"wafy",typeof(wafy),(/"nvar","year","lev","lat","lon"/))
filevardef(ncdf,"wafz",typeof(wafz),(/"nvar","year","lev","lat","lon"/))
else
filevardef(ncdf,"wafx",typeof(wafx),(/"nvar","lev","lat","lon"/))
filevardef(ncdf,"wafy",typeof(wafy),(/"nvar","lev","lat","lon"/))
filevardef(ncdf,"wafz",typeof(wafz),(/"nvar","lev","lat","lon"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"wafx",wafx)
filevarattdef(ncdf,"wafy",wafy)
filevarattdef(ncdf,"wafz",wafz)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->wafx = (/wafx/)
ncdf->wafy = (/wafy/)
ncdf->wafz = (/wafz/)
end

