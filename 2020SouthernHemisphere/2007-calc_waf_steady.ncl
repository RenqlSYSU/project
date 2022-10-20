;***************************************************************************
;calc yearly steady wave activity flux
;                                            by Ql_Ren
;                                           2020/11/10
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1  = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/2020SouthernHemisphere/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","air","hgt"/)
fileout  = path2 + "ERA-Interim_steady_Pb85-year_month.nc"
nmonth   = 3

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
lats =  -5 
latn =  -90
lonl =  0  
lonr = 360
;lev  = (/850,500,200,100/)
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
nlev  = dimsizes(lev)
edd_sf= new((/nyear,nmonth,nlev,nlat,nlon/),float)
wafx  = new((/nyear,nmonth,nlev,nlat,nlon/),float)
wafy  = new((/nyear,nmonth,nlev,nlat,nlon/),float)
wafz  = new((/nyear,nmonth,nlev,nlat,nlon/),float)

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
do nm = 0, nmonth-1, 1
month_s = nm + 3
month_e = nm + 3
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
u     := new((/nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/nyear,ntime,nlev,nlat,nlon/),float)
t     := new((/nyear,ntime,nlev,nlat,nlon/),float)
z     := new((/nyear,ntime,nlev,nlat,nlon/),float)

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
u  := dim_avg_n(u,itime)
v  := dim_avg_n(v,itime)
t  := dim_avg_n(t,itime)
z  := dim_avg_n(z,itime)
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

;the Brunt-Vaisala frequency of atmosphere
theta = pot_temp(lev*100, t, ilev-1, False)
delete(t)
hgt = z/9.8 ;convert unit from m2/s2 to m
NN  = brunt_vaisala_atm(theta, hgt, 0, ilev-1)
delete([/theta,hgt/])
NN = NN*abs(NN)
NN = conform(NN, wgt_areaave(NN, a*cos(lat)*(lon(1)-lon(0)), 1.0, 0), (/0,1/))
;NN = where(abs(NN).lt.0.0000000001, 0.000001, NN)

;stream function
f0 = conform(z,2.0*(2.0*pi/24.0/3600.0)*sin(lat),ilat-1)
;sf = u
;vp = v
;uv2sfvpf ( u, v, sf, vp )
;delete(vp)
sf = z/f0
NN = f0*f0/NN
delete([/z,f0/])

if (any(ismissing(sf))) then
    print("stream function has missing")
end if

delete([/u,v/])
;NN = conform(sf, dim_avg_n(NN,ilon-1), (/0,1,2/))
clm_sf = dim_avg_n(sf,ilon-1)
eddsf = sf - conform(sf, clm_sf, (/0,1,2/))
delete([/clm_sf,sf/])

coslat = conform(eddsf,cos(lat),ilat-1)
;coef  = coslat/2.0
coef   = (coslat*conform(eddsf,lev,ilev-1)/1000.0)/2.0

dsfdx = center_finite_diff_n(eddsf,lon,loncyc,opt,ilon-1)
wafx(:,nm,:,:,:)  = coef*(1.0/a/a/coslat/coslat)*(dsfdx*dsfdx-eddsf*center_finite_diff_n(dsfdx,lon,loncyc,opt,ilon-1))
printMinMax(wafx,False)

dsfdy = center_finite_diff_n(eddsf,lat,cyclic,opt,ilat-1)
wafy(:,nm,:,:,:)  = coef*(1.0/a/a/coslat)*(dsfdx*dsfdy-eddsf*center_finite_diff_n(dsfdy,lon,loncyc,opt,ilon-1))
printMinMax(wafy,False)
delete(dsfdy)

noz    = -sclhgt*log(lev/1000.0)
dsfdz = center_finite_diff_n(eddsf,noz,cyclic,opt,ilev-1)
wafz(:,nm,:,:,:)  = coef*NN*(1.0/a/coslat)*(dsfdx*dsfdz-eddsf*center_finite_diff_n(dsfdz,lon,loncyc,opt,ilon-1))
delete([/coslat,coef,dsfdx,dsfdz/])
printMinMax(wafz,False)
edd_sf(:,nm,:,:,:) = eddsf
end do

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
dimNames = (/"year" ,"month","lev","lat","lon"/);
dimSizes = (/nyear  ,nmonth ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"edd_sf",typeof(edd_sf),(/"year","month","lev","lat","lon"/))
filevardef(ncdf,"wafx",typeof(wafx),(/"year","month","lev","lat","lon"/))
filevardef(ncdf,"wafy",typeof(wafy),(/"year","month","lev","lat","lon"/))
filevardef(ncdf,"wafz",typeof(wafz),(/"year","month","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"edd_sf",edd_sf)
filevarattdef(ncdf,"wafx",wafx)
filevarattdef(ncdf,"wafy",wafy)
filevarattdef(ncdf,"wafz",wafz)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->edd_sf = (/edd_sf/)
ncdf->wafx = (/wafx/)
ncdf->wafy = (/wafy/)
ncdf->wafz = (/wafz/)
end

