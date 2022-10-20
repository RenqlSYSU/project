;***************************************************************************
;regress the anomalous wind, then calc wave activity flux
;                                            by Ql_Ren
;                                           2020/11/10
;******************************************************************************
begin
case = (/"NCEP2","ERA-Interim"/)
nc = 0
if(nc.eq.0) then 
path1  = "/home/ys17-19/data-observation/NCEP_DOE/daily/" ;hgt unit is m 
years = ispan(1979,2018,1)
else
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/" ;hgt unit is m2/s2
years = ispan(1979,2016,1)
end if
path2   = "/home/ys17-19/renql/project/2020SouthernHemisphere/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","air","hgt"/)

;ts = (/-16.33339,-10.07344,-2.129999,-7.827043,26.10205,34.93106,7.237824,-2.88431,-4.954339,-0.844423,-5.536781,-6.130042,-14.13272,4.539642,-12.55697,\
;2.724565,-2.883797,-3.598295,-16.95587,19.35723,13.80637,13.8677,8.82702,-10.5183,-7.378797,-14.55946,-7.588447,5.988917,-8.355371,9.735603,6.808292,\
;-3.840728,7.263536,12.45822,3.979857,-8.741457,-12.33934,2.832654/) ;,2.65007,-2.947284
ts = asciiread(path2+"JET_index_MEJindx_2040N5080E_7918_DJF.txt",-1,"float")
ts := ts(1::)
printVarSummary(ts)
print(ts)

waf_option = 0
if(waf_option.eq.0) then 
fileout  = path2 + case(nc) + "_DJF_steady_TN01-reg.nc"
end if
if(waf_option.eq.1) then 
fileout  = path2 + case(nc) + "_DJF_steady_Pb85-reg2.nc"  ;use streamfunction to calc
end if                                                                                
if(waf_option.eq.2) then                                                              
fileout  = path2 + case(nc) + "_DJF_steady_Pb85_2-reg.nc" ;use anomalous wind to calc
end if
DJF_option = True

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
sclhgt = 7991 ;m 
pi = atan(1.0)*4
print(pi)
years@calendar = "standard"

;area: east Asia
lats =  5  ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/400,300,250, 200,150,100/) 
;lev  = (/ 1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100, 70, 50, 30, 20, 10 /) ;NCEP2
;lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/) ;ERA-Interim

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

if(nc.eq.0) then 
f    = addfile(path1+var_name(0)+"/"+var_name(0)+"."+years(0)+".nc", "r")
else
f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
end if
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
print(nyear)

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
if(DJF_option) then 
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
else
select_time  = ind(time(:,1).ge.3.and.time(:,1).le.5)
end if
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
z     := new((/nyear,ntime,nlev,nlat,nlon/),float)

if(DJF_option) then 
if(nc.eq.0) then 
do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr}))
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr}))
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr}))
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr}))
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
end if
end do
else
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
end if
;--------------------------------------
;DJF
;---------------------------------------------------
u2 := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2 := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
t2 := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
z2 := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    t2(ny,0:30,:,:,:) = t(ny,59:89,:,:,:)
    z2(ny,0:30,:,:,:) = z(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    t2(ny,31:89,:,:,:) = t(ny+1,0:58,:,:,:)
    z2(ny,31:89,:,:,:) = z(ny+1,0:58,:,:,:)
end do
u  := dim_avg_n(u2,itime)
v  := dim_avg_n(v2,itime)
t  := dim_avg_n(t2,itime)
z  := dim_avg_n(z2,itime)

else

if(nc.eq.0) then 
do nt = 0,nyear-1,1
print("read: " + years(nt))
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    t(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    nv = 3
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    z(nt,:,:,:,:)  = short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
end do
else
do nt = 0,nyear-1,1
    print("read: " + years(nt))
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
end do
end if
u  := dim_avg_n(u,itime)
v  := dim_avg_n(v,itime)
t  := dim_avg_n(t,itime)
z  := dim_avg_n(z,itime)
end if
printVarSummary(u)

if (any(ismissing(z))) then
    print("hgt has missing")
end if
if (any(ismissing(u))) then
    print("uwnd has missing")
end if
if (any(ismissing(v))) then
    print("vwnd has missing")
end if
z = linmsg_n(z,-1,ilon-1)
u = linmsg_n(u,-1,ilon-1)
v = linmsg_n(v,-1,ilon-1)
t = linmsg_n(t,-1,ilon-1)

if(waf_option.eq.2) then ;pb85_2

noz    = -sclhgt*log(lev/1000.0)
area_t = wgt_areaave(dim_avg_n_Wrap(t,iyear), a*cos(lat)*(lon(1)-lon(0)), 1.0, 0)
static = center_finite_diff_n(area_t,noz,cyclic,opt,0) + 0.286*area_t/sclhgt
delete([/area_t,noz/])

;edd_u = u - conform(u,dim_avg_n(u,ilon-1),(/0,1,2/))
;edd_v = v - conform(u,dim_avg_n(v,ilon-1),(/0,1,2/))
;edd_t = t - conform(u,dim_avg_n(t,ilon-1),(/0,1,2/))
;edd_z = z - conform(u,dim_avg_n(z,ilon-1),(/0,1,2/))
edd_u = regCoef_n(ts, u, 0, 0) ;narea,neval,nlev,nlat,nlon
edd_v = regCoef_n(ts, v, 0, 0)
edd_t = regCoef_n(ts, t, 0, 0)
edd_z = regCoef_n(ts, z, 0, 0)

tval  := onedtond(edd_z@tval ,(/narea,neval,nlev,nlat,nlon/))
df    := onedtond(edd_z@nptxy,(/narea,neval,nlev,nlat,nlon/))
lon_m := conform(tval,vars&lon,ilon);gitude
printVarSummary(tval)
printVarSummary(df)
b := tval 
b = 0.5
prob  = betainc( df/(df+tval^2), df/2.0, b )
prob  = mask(prob,lon_m.eq.0,False)
delete([/tval,df,lon_m/])

f0     = conform(edd_z, 2.0*(2.0*pi/24.0/3600.0)*sin(lat), ilat)
edd_sf = edd_z/f0
delete([/u,v,t,z,f0/])
static := conform(edd_u, static, (/ilev/))

coef1  = conform(edd_u,2.0*a*(2.0*pi/24.0/3600.0)*sin(2.0*lat),ilat)
coef2  = conform(edd_u,2.0*(2.0*pi/24.0/3600.0)*sin(lat),ilat)/static
coef3  = conform(edd_u,cos(lat),ilat)*conform(edd_u,lev/1000.0,ilev)
delete(static)

wafy = coef3 * (-edd_v*edd_u/2.0 + (1.0/coef1)*center_finite_diff_n(edd_u*edd_z,lon,loncyc,opt,ilon)) 
delete(edd_u)
wafz = coef3 * ( edd_v*edd_t/2.0 - (1.0/coef1)*center_finite_diff_n(edd_t*edd_z,lon,loncyc,opt,ilon))*coef2
delete([/edd_t,coef2/])
wafx = coef3 * ( edd_v*edd_v/2.0 - (1.0/coef1)*center_finite_diff_n(edd_v*edd_z,lon,loncyc,opt,ilon)) 
delete([/edd_v,edd_z,coef1,coef3/])
printMinMax(wafx,False)
printMinMax(wafy,False)
printMinMax(wafz,False)

else

;the Brunt-Vaisala frequency of atmosphere
theta = pot_temp(lev*100, dim_avg_n_Wrap(t,iyear), ilev-2, False)
delete(t)

if(nc.eq.1)
hgt = dim_avg_n_Wrap(z,iyear)/9.8 ;convert unit from m2/s2 to m
else
hgt = dim_avg_n_Wrap(z,iyear) ;convert unit from m2/s2 to m
z   = z*9.8  ;convert unit from m to m2/s2
end if
NN  = brunt_vaisala_atm(theta, hgt, 0, ilev-2)
delete([/theta,hgt/])
NN := wgt_areaave(NN*abs(NN), a*cos(lat)*(lon(1)-lon(0)), 1.0, 0) ;(lev)
;NN = where(abs(NN).lt.0.0000000001, 0.000001, NN)

;stream function
;sf = u
;vp = v
;uv2sfvpf ( u, v, sf, vp )
;delete(vp)

;z = z - conform(z, dim_avg_n_Wrap(z,ilon-1), (/0,1,2/))
edd_sf = regCoef_n(ts, z, 0, 0) ;narea,neval,nlev,nlat,nlon
printVarSummary(edd_sf)
;edd_sf = edd_sf - conform(edd_sf, dim_avg_n_Wrap(edd_sf,ilon), (/0,1,2,3/))
delete(z)
tval  := onedtond(edd_sf@tval ,(/nlev,nlat,nlon/))
df    := onedtond(edd_sf@nptxy,(/nlev,nlat,nlon/))
lon_m := conform(tval,vars&lon,ilon-2);gitude
printVarSummary(tval)
printVarSummary(df)
b := tval 
b = 0.5
prob  = betainc( df/(df+tval^2), df/2.0, b )
prob  = mask(prob,lon_m.eq.0,False)
delete([/tval,df,lon_m/])

f0     = conform(edd_sf, 2.0*(2.0*pi/24.0/3600.0)*sin(lat), ilat-2)
edd_sf = edd_sf/f0

NN := conform(edd_sf, NN, (/ilev-2/))
NN = f0*f0/NN
delete(f0)

if (any(ismissing(edd_sf))) then
    print("stream function has missing")
end if

if(waf_option.eq.0) then ;TN01
print("")
u  := dim_avg_n(u ,iyear)
v  := dim_avg_n(v ,iyear)

dsfdx = center_finite_diff_n(edd_sf,lon,loncyc,opt,ilon-2)
dsfdy = center_finite_diff_n(edd_sf,lat,cyclic,opt,ilat-2)

coslat = conform(edd_sf,cos(lat),ilat-2)
;coef   = coslat/2.0/sqrt(u^2 + v^2)
coef   = (coslat*conform(u,lev/1000.0,ilev-2))/2.0/sqrt(u^2 + v^2)
wafx  = coef*((u/a/a/coslat/coslat)*(dsfdx*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lon,loncyc,opt,ilon-2)) \
             +(v/a/a/coslat       )*(dsfdy*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lat,cyclic,opt,ilat-2)))
printMinMax(wafx,False)

wafy  = coef*((u/a/a/coslat)*(dsfdx*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lon,loncyc,opt,ilon-2)) \
             +(v/a/a       )*(dsfdy*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lat,cyclic,opt,ilat-2)))
printMinMax(wafy,False)

noz   = -sclhgt*log(lev/1000.0)
dsfdz = center_finite_diff_n(edd_sf,noz,cyclic,opt,ilev-2)
;dsfdz = center_finite_diff_n(edd_sf,lev*100,cyclic,opt,ilev-2)
wafz  = coef*NN*((u/a/coslat)*(dsfdx*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lon,loncyc,opt,ilon-2)) \
                +(v/a       )*(dsfdy*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lat,cyclic,opt,ilat-2)))
delete([/coslat,coef,dsfdx,dsfdy,dsfdz,u,v/])
printMinMax(wafz,False)
end if

if(waf_option.eq.1) then
delete([/u,v/])
;clm_sf = dim_avg_n(sf,ilon-1)
;edd_sf = sf - conform(sf, clm_sf, (/0,1,2/))
;delete([/clm_sf,sf/])
;NN = conform(edd_sf, dim_avg_n(NN,ilon-1), (/0,1,2/))

coslat = conform(edd_sf,cos(lat),ilat)
;coef   = coslat/2.0
coef   = (coslat*conform(edd_sf,lev/1000.0,ilev))/2.0

dsfdx = center_finite_diff_n(edd_sf,lon,loncyc,opt,ilon)
wafx  = coef*(1.0/a/a/coslat/coslat)*(dsfdx*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lon,loncyc,opt,ilon))
printMinMax(wafx,False)

dsfdy = center_finite_diff_n(edd_sf,lat,cyclic,opt,ilat)
wafy  = coef*(1.0/a/a/coslat)*(dsfdx*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lon,loncyc,opt,ilon))
printMinMax(wafy,False)
delete(dsfdy)

noz    = -sclhgt*log(lev/1000.0)
dsfdz = center_finite_diff_n(edd_sf,noz,cyclic,opt,ilev)
wafz  = coef*NN*(1.0/a/coslat)*(dsfdx*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lon,loncyc,opt,ilon))
delete([/coslat,coef,dsfdx,dsfdz/])
printMinMax(wafz,False)
end if
end if

edd_sf@long_name = "a perturbation streamfunction based on winter climatology"
edd_sf@units     = "m2/s" 
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
dimNames = (/"lev","lat","lon"/);
dimSizes = (/nlev ,nlat ,nlon/) ;
dimUnlim = (/False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"edd_sf",typeof(edd_sf),(/"lev","lat","lon"/))
filevardef(ncdf,"prob"  ,typeof(prob)  ,(/"lev","lat","lon"/))
filevardef(ncdf,"wafx",typeof(wafx),(/"lev","lat","lon"/))
filevardef(ncdf,"wafy",typeof(wafy),(/"lev","lat","lon"/))
filevardef(ncdf,"wafz",typeof(wafz),(/"lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"edd_sf",edd_sf)
filevarattdef(ncdf,"prob"  ,prob)
filevarattdef(ncdf,"wafx",wafx)
filevarattdef(ncdf,"wafy",wafy)
filevarattdef(ncdf,"wafz",wafz)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->edd_sf = (/edd_sf/)
ncdf->prob   = (/prob/)
ncdf->wafx = (/wafx/)
ncdf->wafy = (/wafy/)
ncdf->wafz = (/wafz/)
end

