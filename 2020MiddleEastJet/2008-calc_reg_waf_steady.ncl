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
waf_option = 1
if(waf_option.eq.0) then 
fileout  = path2 + "ERA-Interim_DJF_steady_TN01-reg.nc"
end if
if(waf_option.eq.1) then 
fileout  = path2 + "ERA-Interim_DJF_steady_Pb85-reg3.nc"
end if
if(waf_option.eq.2) then 
fileout  = path2 + "ERA-Interim_DJF_steady_Pb85_2-reg.nc"
end if

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
sclhgt = 7991 ;m 
pi = atan(1.0)*4
print(pi)
years = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats =  5  ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)
;lev  = (/925, 200,175,150,125,100/)

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

;------------------EOF time series------------------------
neval = 3
area  = (/"ME","EA","NA"/);,"GB"
narea = dimsizes(area)
ts    = new((/narea,neval,nyear-1/),float)
eval  = new((/narea,neval/),float)
pcvar = new((/narea,neval/),float)
sig   = new((/narea,neval/),logical)
do na = 0, narea-1,1
    f  = addfile(path2 + "eof_u_6kinds-" + area(na) + ".nc","r")
    eof0 := f->eof0
    eval(na,:)  = eof0@eval(0:(neval-1))
    pcvar(na,:) = eof0@pcvar(0:(neval-1))
    ts(na,:,:)  = f->eof_ts(0, 0:(neval-1), 0:(nyear-2))
    sig(na,:)   = eofunc_north(pcvar(na,:),(nyear-1),False)
end do
ts = - ts/conform(ts,sqrt(eval),(/0,1/))
delete(eof0)
;--------------------------------------------------------

if(waf_option.eq.2) then

noz    = -sclhgt*log(lev/1000.0)
area_t = wgt_areaave(dim_avg_n_Wrap(t,iyear), a*cos(lat)*(lon(1)-lon(0)), 1.0, 0)
static = center_finite_diff_n(area_t,noz,cyclic,opt,ilev-1) + 0.286*area_t/sclhgt
delete([/area_t,noz/])

;edd_u = u - conform(u,dim_avg_n(u,ilon-1),(/0,1,2/))
;edd_v = v - conform(u,dim_avg_n(v,ilon-1),(/0,1,2/))
;edd_t = t - conform(u,dim_avg_n(t,ilon-1),(/0,1,2/))
;edd_z = z - conform(u,dim_avg_n(z,ilon-1),(/0,1,2/))
edd_u = regCoef_n(ts, u, 2, 0) ;narea,neval,nlev,nlat,nlon
edd_v = regCoef_n(ts, v, 2, 0)
edd_t = regCoef_n(ts, t, 2, 0)
edd_z = regCoef_n(ts, z, 2, 0)
delete([/u,v,t,z/])
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
hgt = dim_avg_n_Wrap(z,iyear)/9.8 ;convert unit from m2/s2 to m
NN  = brunt_vaisala_atm(theta, hgt, 0, ilev-2)
delete([/theta,hgt/])
NN := wgt_areaave(NN*abs(NN), a*cos(lat)*(lon(1)-lon(0)), 1.0, 0) ;(lev)
;NN = where(abs(NN).lt.0.0000000001, 0.000001, NN)

;stream function
;sf = u
;vp = v
;uv2sfvpf ( u, v, sf, vp )
;delete(vp)

z = z - conform(z, dim_avg_n_Wrap(z,ilon-1), (/0,1,2/))
edd_sf = regCoef_n(ts, z, 2, 0) ;narea,neval,nlev,nlat,nlon
printVarSummary(edd_sf)
delete(z)
tval  := onedtond(edd_sf@tval ,(/narea,neval,nlev,nlat,nlon/))
df    := onedtond(edd_sf@nptxy,(/narea,neval,nlev,nlat,nlon/))
lon_m := conform(tval,vars&lon,ilon);gitude
printVarSummary(tval)
printVarSummary(df)
b := tval 
b = 0.5
prob  = betainc( df/(df+tval^2), df/2.0, b )
prob  = mask(prob,lon_m.eq.0,False)
delete([/tval,df,lon_m/])

f0     = conform(edd_sf, 2.0*(2.0*pi/24.0/3600.0)*sin(lat), ilat)
edd_sf = edd_sf/f0

NN := conform(edd_sf, NN, (/ilev/))
NN = f0*f0/NN
delete(f0)

if (any(ismissing(edd_sf))) then
    print("stream function has missing")
end if

;==========================================================
;calc forcing factors: Qeddy, Qd, A, dimension is 5d
;============================================================
if(waf_option.eq.0) then
print("")
u  = conform(edd_sf, dim_avg_n(u ,iyear), (/ilev,ilat,ilon/))
v  = conform(edd_sf, dim_avg_n(v ,iyear), (/ilev,ilat,ilon/))
;NN = conform(edd_sf, dim_avg_n(NN,iyear), (/ilev,ilat,ilon/))
;clm_sf = dim_avg_n(sf,iyear)
;edd_sf = sf - conform(sf, clm_sf, (/1,2,3/))
;delete([/clm_sf/])

dsfdx = center_finite_diff_n(edd_sf,lon,loncyc,opt,ilon)
dsfdy = center_finite_diff_n(edd_sf,lat,cyclic,opt,ilat)

coslat = conform(edd_sf,cos(lat),ilat)
coef   = coslat/2.0/sqrt(u^2 + v^2)
;coef   = (coslat*conform(u,lev,ilev)/1000.0)/2.0/sqrt(u^2 + v^2)
wafx  = coef*((u/a/a/coslat/coslat)*(dsfdx*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lon,loncyc,opt,ilon)) \
             +(v/a/a/coslat       )*(dsfdy*dsfdx-edd_sf*center_finite_diff_n(dsfdx,lat,cyclic,opt,ilat)))
printMinMax(wafx,False)

wafy  = coef*((u/a/a/coslat)*(dsfdx*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lon,loncyc,opt,ilon)) \
             +(v/a/a       )*(dsfdy*dsfdy-edd_sf*center_finite_diff_n(dsfdy,lat,cyclic,opt,ilat)))
printMinMax(wafy,False)

dsfdz = center_finite_diff_n(edd_sf,lev*100,cyclic,opt,ilev)
wafz  = coef*NN*((u/a/coslat)*(dsfdx*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lon,loncyc,opt,ilon)) \
                +(v/a       )*(dsfdy*dsfdz-edd_sf*center_finite_diff_n(dsfdz,lat,cyclic,opt,ilat)))
delete([/coslat,coef,dsfdx,dsfdy,dsfdz/])
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
coef   = (coslat*conform(edd_sf,lev,ilev)/1000.0)/2.0

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

edd_sf@long_name = "a perturbation streamfunction based on winter climatology"
edd_sf@units     = "m2/s" 
end if

wafx@long_name = "x wave activity flux"
wafx@units   = "m2/s2" 
wafx@pcvar   = ndtooned(pcvar)
wafx@eval    = ndtooned(eval)
wafx@sig     = ndtooned(sig)
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
dimNames = (/"area","neval","lev","lat","lon"/);
dimSizes = (/narea ,neval  ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(waf_option.le.1) then 
filevardef(ncdf,"edd_sf",typeof(edd_sf),(/"area","neval","lev","lat","lon"/))
filevardef(ncdf,"prob"  ,typeof(prob)  ,(/"area","neval","lev","lat","lon"/))
end if
filevardef(ncdf,"wafx",typeof(wafx),(/"area","neval","lev","lat","lon"/))
filevardef(ncdf,"wafy",typeof(wafy),(/"area","neval","lev","lat","lon"/))
filevardef(ncdf,"wafz",typeof(wafz),(/"area","neval","lev","lat","lon"/))
filevardef(ncdf,"eval" ,typeof(eval) ,(/"area","neval"/))
filevardef(ncdf,"pcvar",typeof(pcvar),(/"area","neval"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
if(waf_option.le.1) then 
filevarattdef(ncdf,"edd_sf",edd_sf)
filevarattdef(ncdf,"prob"  ,prob)
end if
filevarattdef(ncdf,"wafx",wafx)
filevarattdef(ncdf,"wafy",wafy)
filevarattdef(ncdf,"wafz",wafz)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
if(waf_option.le.1) then 
ncdf->edd_sf = (/edd_sf/)
ncdf->prob   = (/prob/)
end if
ncdf->wafx = (/wafx/)
ncdf->wafy = (/wafy/)
ncdf->wafz = (/wafz/)
ncdf->eval  = (/eval/)
ncdf->pcvar = (/pcvar/)
end

