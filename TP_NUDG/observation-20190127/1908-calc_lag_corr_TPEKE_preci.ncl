;******************************************************************************
;1. calculate TP area average vertically integrated EKE
;2. calculate the correlation between
;                                            by Ql_Ren
;                                           2019/08/19
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc"
path5   = "/home/ys17-19/data-observation/ERA-interim/surface/rain_daily/rain.daily.interim." ;+years(nt)+".nc"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc"
var_name = (/"uwnd","vwnd","r"/)
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter

;parameter for bw_bandpass_filter
ca = 7.0 ;day
cb = 2.5 

if(option.eq.0) then
    ;fileout  = path2 + "mdata/ERA-Interim_corr_intEKE_preci.nc"
    fileout  = path2 + "mdata/ERA-Interim_lag_corr_DJFintEKE_JJApreci.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_lag_corr_DJFintEKE-" + cb + "-" + ca + "day_JJApreci.nc"
    ;fileout  = path2 + "mdata/ERA-Interim_corr_intEKE-" + cb + "-" + ca + "day_preci.nc"
end if
g  = 9.8 ;m/(s*s)
years    = ispan(1979,2005,1)
years@calendar = "standard"

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"Tibet"/)
lats = (/25/)
latn = (/40/)
lonl = (/70/)
lonr = (/105/)
lev  = (/1000,925,850,700,600,500,400,300,200/)

;area of field variables,larger then the plan
flats = -90 ;-15
flatn =  90 ;55
flonl =  0  ;30
flonr = 360 ;150

f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
printVarSummary(time)

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats:latn},{lonl:lonr}))
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f     = addfile(path5 + years(0) + ".nc","r")
vars1 = f->$var_name(2)$(0,{lats:latn},{lonl:lonr})
ntime = 365
nyear = dimsizes(years)
nlat = dimsizes(vars1&lat)
nlon = dimsizes(vars1&lon)
nlev = dimsizes(lev)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)

vars  = f->$var_name(2)$(0,{flats:flatn},{flonl:flonr})
nflat = dimsizes(vars&lat)
nflon = dimsizes(vars&lon)
EKEts = new((/nyear,ntime/),float)
fvars = new((/nyear,ntime,nflat,nflon/),float)
corr  = new((/nflat,nflon/),float)

;==========================================================
;calculate the EKE index and read the field variable
;============================================================
index  = new((/ntime/),integer)
index1 = new((/ntime/),integer)
index(0:58)   = ispan(0,58,1)
index(59:364) = ispan(60,365,1)

f_sfcp  = addfile(path4,"r")
start = 0
do nt = 0,nyear-1,1
print("read: " + years(nt) + ", start = " + start)
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path5+years(nt)+".nc","r")
    fvars(nt,:,:,:)  = f->$var_name(nv)$(index,{flats:flatn},{flonl:flonr})
    space = 366
    index1(0:58)   = ispan(start,(start+58),1)
    index1(59:364) = ispan((start+60),(start+space-1),1)
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path5+years(nt)+".nc","r")
    fvars(nt,:,:,:)  = f->$var_name(nv)$(:,{flats:flatn},{flonl:flonr})
    space = 365
    ps(nt,:,:,:)   = f_sfcp->sp(start:(start+space-1):1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
fvars = where(fvars.lt.0,0,fvars)
printVarSummary(ps)
printVarSummary(u)
printVarSummary(dp)

if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb
    opt = False
    filt_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    filt_v = bw_bandpass_filter(v,fca,fcb,opt,1)
    var   = dim_sum_n((filt_v*filt_v + filt_u*filt_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
else
    ste_u = conform(u,dim_avg_n(u,0),(/1,2,3,4/))
    ste_v = conform(v,dim_avg_n(v,0),(/1,2,3,4/))
    edd_u = u - conform(u,dim_avg_n(u,0),(/1,2,3,4/))
    edd_v = v - conform(v,dim_avg_n(v,0),(/1,2,3,4/))
    var   = dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
end if
var = mask(var,phis.gt.1500,True) 
EKEts = wgt_areaave_Wrap(var,1.0,1.0,0)
printVarSummary(EKEts)

select_time1 = ind(time(:,1).eq.12.or.time(:,1).le.2)
select_time2 = ind(time(:,1).ge.6.and.time(:,1).le.8)
season_corr  = escorc_n(dim_avg_n_Wrap(EKEts(:,select_time1),1),dim_avg_n_Wrap(fvars(:,select_time2,:,:),1),0,0)
printVarSummary(fvars)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the correlation between Asia rainfall and TP area average EKE "
fileattdef(ncdf,fileAtt)  

;define all coordinate
dimNames = (/"season","lat","lon"/);
dimSizes = (/  1     ,nflat,nflon/) ;
dimUnlim = (/False   ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat" ,typeof(fvars&lat),getvardims(fvars&lat))
filevardef(ncdf,"lon" ,typeof(fvars&lon),getvardims(fvars&lon))
filevardef(ncdf,"season_corr" ,typeof(season_corr) ,(/"season","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",fvars&lat)
filevarattdef(ncdf,"lon",fvars&lon)
filevarattdef(ncdf,"season_corr",season_corr)

ncdf->season = (/1212/)
ncdf->lat  = (/fvars&lat/)
ncdf->lon  = (/fvars&lon/)
ncdf->season_corr  = (/season_corr/)
end

