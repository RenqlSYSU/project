;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc"
var_name = (/"uwnd","vwnd"/)

year_option = True;False ;
fileout  = path2 + "ERA-Interim_DJF_intEKE-year.nc"

g  = 9.8 ;m/(s*s)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats =  10 ;-20
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
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
ntime = dimsizes(select_time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)

;==========================================================
;read the data
;============================================================
f_sfcp  = addfile(path4,"r")
start = 0
do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    space = 366
    index1 = select_time1 + start
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    space = 365
    index1 = select_time + start
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do

;--------------------------------------
;DJF
;---------------------------------------------------
u2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
dp2    := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:) = u(ny,59:89,:,:,:)
    v2(ny,0:30,:,:,:) = v(ny,59:89,:,:,:)
    dp2(ny,0:30,:,:,:)=dp(ny,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
    dp2(ny,31:89,:,:,:)=dp(ny+1,0:58,:,:,:)
end do
u := u2
v := v2
dp:= dp2
printVarSummary(u)

;--------------------------------------
;Calculate seasonal averages and anomalies based on seasonal averages
;---------------------------------------------------
itime = 1
ca = 10.0
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
delete([/cf,term,fcf,fca,hgh_cf,low_cf/])

clm_u = dim_avg_n(u,1)
clm_v = dim_avg_n(v,1) 
edd_u = u - conform(u, clm_u, (/0,2,3,4/))
edd_v = v - conform(v, clm_v, (/0,2,3,4/))   
printVarSummary(clm_u)
delete([/u,v/])
delete([/clm_u,clm_v/])

var = dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
EKE = dim_avg_n(var,1) 
EKE@long_name = "DJF Integrated transient eddies kinetic energy"
EKE@units     = "J" 

var    = dim_sum_n((hgh_v*hgh_v + hgh_u*hgh_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
hghEKE = dim_avg_n(var,1) 
hghEKE@long_name = "DJF Integrated high frequency transient eddies kinetic energy"
hghEKE@units     = "J" 

var    = dim_sum_n((low_v*low_v + low_u*low_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
lowEKE = dim_avg_n(var,1) 
lowEKE@long_name = "DJF Integrated low frequecy transient eddies kinetic energy"
lowEKE@units     = "J" 
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
fileatt@discription   = "DJF" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year" ,"lat","lon"/);
dimSizes = (/nyear-1,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"intEKE",typeof(EKE),(/"year","lat","lon"/))
filevardef(ncdf,"hghEKE",typeof(hghEKE),(/"year","lat","lon"/))
filevardef(ncdf,"lowEKE",typeof(lowEKE),(/"year","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"intEKE",EKE)
filevarattdef(ncdf,"hghEKE",hghEKE)
filevarattdef(ncdf,"lowEKE",lowEKE)

ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->intEKE = (/EKE/)
ncdf->hghEKE = (/hghEKE/)
ncdf->lowEKE = (/lowEKE/)
end

