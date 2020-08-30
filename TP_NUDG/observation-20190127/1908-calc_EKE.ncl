;***************************************************************************
;use the post-processing data to calculate climatology daily transient eddies
;1. do bandpass filtering;  2. calc transient eddies kinetic energy (EkE)
;3. calc climatology daily EKE
;4. store EKE on each level
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/01/25
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc"
var_name = (/"uwnd","vwnd"/)
option   = 1 ;0 means remove annual cycle, 1 means bandpass filter

;parameter for bw_bandpass_filter
ca = 90.0 ;day
cb = 10.0 

if(option.eq.0) then
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_EKE.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_EKE-" + cb + "-" + ca + "day.nc"
end if

g  = 9.8 ;m/(s*s)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats =  10 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600,500,400,300,200,150,100,50,30,10,5/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
    t   = f_model->U(:,0,{lats:latn},{lonl:lonr})
printVarSummary(t)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,0,{lats:latn},{lonl:lonr})
ntime = 365
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)
var_month  = new((/nyear,12 ,nlev,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlev,nlat,nlon/),float)
month_ave  = new((/      12 ,nlev,nlat,nlon/),float)
season_ave = new((/       4 ,nlev,nlat,nlon/),float)

;==========================================================
;read the data
;============================================================
index  = new((/ntime/),integer)
index1 = new((/ntime/),integer)
index(0:58)   = ispan(0,58,1)
index(59:364) = ispan(60,365,1)

do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
end if
end do
printVarSummary(u)

if(option.eq.1) then
fca = 1.0/ca
fcb = 1.0/cb
opt = False
filt_u = bw_bandpass_filter(u,fca,fcb,opt,1)
filt_v = bw_bandpass_filter(v,fca,fcb,opt,1)
print(filt_u(0,0:10,0,30,30))
var   = (filt_v*filt_v + filt_u*filt_u)/2.0   ;get (nyear,ntime,nlat,nlon)
else
ste_u = conform(u,dim_avg_n(u,0),(/1,2,3,4/))
ste_v = conform(v,dim_avg_n(v,0),(/1,2,3,4/))
edd_u = u - conform(u,dim_avg_n(u,0),(/1,2,3,4/))
edd_v = v - conform(v,dim_avg_n(v,0),(/1,2,3,4/))
var   = (edd_v*edd_v + edd_u*edd_u)/2.0   ;get (nyear,ntime,nlat,nlon)
end if
printVarSummary(var)
print(var(0,1,8,50,20:30))

EKE   = dim_avg_n(var,0)  ;daily 

EKE@long_name = "Daily Integrated transient eddies kinetic energy"
;EKE@long_name = "Integrated all kinetic energy"
;EKE@long_name = "Integrated steady kinetic energy"
EKE@units     = "J" 

time = cd_calendar(f_model->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

;----------calculate the season average-----------------------
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    var_season(:,ns,:,:,:) = ( var_month(:,11  ,:,:,:)*mdays(11)   + \
                               var_month(:,nm  ,:,:,:)*mdays(nm)   + \
                               var_month(:,nm+1,:,:,:)*mdays(nm+1))/(mdays(11) + mdays(nm)+ mdays(nm+1))
    else
    var_season(:,ns,:,:,:) = ( var_month(:,nm-1,:,:,:)*mdays(nm-1) + \
                               var_month(:,nm  ,:,:,:)*mdays(nm)   + \
                               var_month(:,nm+1,:,:,:)*mdays(nm+1))/(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
    end if
end do
month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
printVarSummary(month_ave)
printVarSummary(season_ave)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
;fileatt@discription   = "the integration is based on 1000,925,850,700,600,500hPa "
fileatt@discription   = "the integration is based on 1000,925,850,700,600,500,400,300,200hPa "
;fileatt@discription   = "the integration is based on 500,400,300,200hPa "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","season","time","lev","lat","lon"/);
dimSizes = (/  12   ,  4     ,ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"time",typeof(t&time),getvardims(t&time))
filevardef(ncdf,"lev" ,typeof(u&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
;filevardef(ncdf,"SKE",typeof(EKE),(/"time","lat","lon"/))
filevardef(ncdf,"EKE",typeof(EKE),(/"time","lev","lat","lon"/))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",t&time)
filevarattdef(ncdf,"lev",u&level)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"EKE",EKE)
;filevarattdef(ncdf,"SKE",EKE)
filevarattdef(ncdf,"month_ave",EKE)
filevarattdef(ncdf,"season_ave",EKE)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->time = (/t&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->EKE  = (/EKE/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
end

