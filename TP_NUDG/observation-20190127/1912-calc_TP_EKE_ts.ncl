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
var_name = (/"uwnd","vwnd"/)
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
na = 0

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

g  = 9.8 ;m/(s*s)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area of index, India,Tibet,Philippines,Maritime_continent
area = (/"TP"/)
lats = (/25 /)
latn = (/40 /)
lonl = (/70 /)
lonr = (/105/)
lev  = (/1000,975,950,925,900, 875,850,825,800,775, 750,700,650,600,550, 500,450,400,350,300, 250,225,200,175,150, 125,100/)

if(option.eq.0) then
    fileout  = path2 + "mdata/ERA-Interim_EKE_" + area(na) + "ts.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_EKE_" + area(na) + "ts-" + cb + "-" + ca + "day.nc"
end if

;==================================================================================
;read data
;=============================================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
printVarSummary(time)

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)}))
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f     = addfile(path5 + years(0) + ".nc","r")
vars  = f->$var_name(2)$(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
ntime = 365 
nyear = dimsizes(years)
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nlev = dimsizes(lev)
u    = new((/nyear,ntime,nlev,nlat,nlon/),float)
v    = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp   = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps   = new((/nyear,ntime,nlat,nlon/),float)
ts   = new((/nyear,ntime/),float)
month_corr     = new((/12,nflat,nflon/),float)
season_corr    = new((/4 ,nflat,nflon/),float)

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

seasons  = (/"DJF","MAM","JJA","SON"/)
months   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    term = dim_avg_n_Wrap(EKEts(:,select_time),1)
    month_corr(nm,:,:) = escorc_n(term,dim_avg_n_Wrap(fvars(:,select_time,:,:),1),0,0)
    print(months(nm))
    print(term)
end do

;mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
;mdays    = (/90,92,92,91/) ;DJF,MAM,JJA,SON
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
        select_time := ind((time(:,1).eq.12).or.(time(:,1).le.2))
    else
        select_time := ind((time(:,1).le.nm).or.(time(:,1).ge.(nm+2)))
    end if
    term = dim_avg_n_Wrap(EKEts(:,select_time),1)
    season_corr(ns,:,:) = escorc_n(term,dim_avg_n_Wrap(fvars(:,select_time,:,:),1),0,0)
    print(seasons(ns))
    print(term)
end do
printVarSummary(fvars)
printVarSummary(month_corr)
printVarSummary(season_corr)

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
dimNames = (/"month","season","lat","lon"/);
dimSizes = (/  12   ,  4     ,nflat,nflon/) ;
dimUnlim = (/False  ,False   ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat" ,typeof(fvars&lat),getvardims(fvars&lat))
filevardef(ncdf,"lon" ,typeof(fvars&lon),getvardims(fvars&lon))
filevardef(ncdf,"month_corr"  ,typeof(month_corr)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_corr" ,typeof(season_corr) ,(/"season","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",fvars&lat)
filevarattdef(ncdf,"lon",fvars&lon)
filevarattdef(ncdf,"month_corr" ,month_corr)
filevarattdef(ncdf,"season_corr",season_corr)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->lat  = (/fvars&lat/)
ncdf->lon  = (/fvars&lon/)
ncdf->month_corr   = (/month_corr/)
ncdf->season_corr  = (/season_corr/)
end

