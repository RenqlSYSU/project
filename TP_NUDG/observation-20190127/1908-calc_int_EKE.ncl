;***************************************************************************
;use the post-processing data to calculate climatology daily transient eddies
;1. remove annual cycle;  2. calc transient eddies kinetic energy (EkE)
;3. calc climatology daily EKE
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
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
filter   = True ;use filtered daily climatology

;parameter for bw_bandpass_filter
ca = 90.0 ;day
cb = 2.5

if(option.eq.0) then
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_uppEKE_filter.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_uppEKE3-" + cb + "-" + ca + "day.nc"
end if
g  = 9.8 ;m/(s*s)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats =  0 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600,\
         500 ,400,350,300,250,\
         200 /) ;,150,100hPa

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
var_month  = new((/nyear,12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)

;==========================================================
;read the data
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
    space = 365
    ps(nt,:,:,:)   = f_sfcp->sp(start:(start+space-1):1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
printVarSummary(ps)
printVarSummary(u)
printVarSummary(dp)

if(option.eq.1) then
;u := reshape(u,(/nyear*ntime,nlev,nlat,nlon/))
;v := reshape(v,(/nyear*ntime,nlev,nlat,nlon/))
fca = 1.0/ca
fcb = 1.0/cb
opt = False
filt_u = bw_bandpass_filter(u,fca,fcb,opt,1)
filt_v = bw_bandpass_filter(v,fca,fcb,opt,1)
delete([/u,v/])
;filt_u := reshape(filt_u,(/nyear,ntime,nlev,nlat,nlon/))
;filt_v := reshape(filt_v,(/nyear,ntime,nlev,nlat,nlon/))

;u := reshape(u,(/nyear*ntime,nlev,nlat,nlon/))
;v := reshape(v,(/nyear*ntime,nlev,nlat,nlon/))
;fca = 1.0/ca
;fcb = 1.0/cb
;fcf = fspan(0, 0.5, ntime*nyear/2) 
;filt_cf = ind(fcf.le.fca.or.fcf.ge.fcb)
;cf = ezfftf_n(u,0)
;cf(:,filt_cf,:,:,:) = 0.0  
;filt_u = ezfftb_n(cf,0,1)
;cf = ezfftf_n(v,0)
;cf(:,filt_cf,:,:,:) = 0.0  
;filt_v = ezfftb_n(cf,0,1)
;filt_u := reshape(filt_u,(/nyear,ntime,nlev,nlat,nlon/))
;filt_v := reshape(filt_v,(/nyear,ntime,nlev,nlat,nlon/))
;
;fca = 1.0/ca
;fcb = 1.0/cb
;fcf = fspan(0, 0.5, ntime/2) 
;filt_cf = ind(fcf.le.fca.or.fcf.ge.fcb)
;cf = ezfftf_n(u,1)
;cf(:,:,filt_cf,:,:,:) = 0.0  
;filt_u = ezfftb_n(cf,0,2)
;cf = ezfftf_n(v,1)
;cf(:,:,filt_cf,:,:,:) = 0.0  
;filt_v = ezfftb_n(cf,0,2)

;var   = dim_sum_n((filt_v*filt_v + filt_u*filt_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
var   = dim_sum_n((filt_v(:,:,5:10,:,:)*filt_v(:,:,5:10,:,:) + filt_u(:,:,5:10,:,:)*filt_u(:,:,5:10,:,:))*dp(:,:,5:10,:,:),2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
delete([/filt_u,filt_v,dp/])

else

if(filter) then 
cf = ezfftf_n(dim_avg_n(u,0),0)
cf(:,4:,:,:,:)  = 0
ste_u = ezfftb_n(cf,cf@xbar,1)
printVarSummary(cf)
cf = ezfftf_n(dim_avg_n(v,0),0)
cf(:,4:,:,:,:)  = 0
ste_v = ezfftb_n(cf,cf@xbar,1)
printVarSummary(ste_v)
delete(cf)
edd_u = u - conform(u,ste_u,(/1,2,3,4/))
edd_v = v - conform(v,ste_v,(/1,2,3,4/))
delete([/ste_u,ste_v/])
else
;ste_u = conform(u,dim_avg_n(u,0),(/1,2,3,4/))
;ste_v = conform(v,dim_avg_n(v,0),(/1,2,3,4/))
edd_u = u - conform(u,dim_avg_n(u,0),(/1,2,3,4/))
edd_v = v - conform(v,dim_avg_n(v,0),(/1,2,3,4/))
end if
delete([/u,v/])
;var   = dim_sum_n((ste_v*ste_v + ste_u*ste_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
;var   = dim_sum_n((v*v + u*u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
;var   = dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
var   = dim_sum_n((edd_v(:,:,5:10,:,:)*edd_v(:,:,5:10,:,:) + edd_u(:,:,5:10,:,:)*edd_u(:,:,5:10,:,:))*dp(:,:,5:10,:,:),2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
delete([/edd_u,edd_v,dp/])
end if

EKE   = dim_avg_n(var,0) 
;EKE_vari = dim_variance_n(var,0) 

EKE@long_name = "Daily Integrated transient eddies kinetic energy"
;EKE@long_name = "Integrated all kinetic energy"
;EKE@long_name = "Integrated steady kinetic energy"
EKE@units     = "J" 

time = cd_calendar(f_model->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

;----------calculate the season average-----------------------
mdays    = (/31,28,31,30,31,30,31,31,30,31,30,31/)
do ns = 0,3,1
    nm = ns*3 ;nm = 0,3,6,9
    if(nm.eq.0) then
    var_season(:,ns,:,:) = ( var_month(:,11  ,:,:)*mdays(11)   + \
                             var_month(:,nm  ,:,:)*mdays(nm)   + \
                             var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(11) + mdays(nm)+ mdays(nm+1))
    else
    var_season(:,ns,:,:) = ( var_month(:,nm-1,:,:)*mdays(nm-1) + \
                             var_month(:,nm  ,:,:)*mdays(nm)   + \
                             var_month(:,nm+1,:,:)*mdays(nm+1))/(mdays(nm-1) + mdays(nm)+ mdays(nm+1))
    end if
end do
month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
season_vari = dim_variance_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
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
fileatt@discription   = "rmvannual cycle"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","season","time","lat","lon"/);
dimSizes = (/  12   ,  4     ,ntime ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
;filevardef(ncdf,"time",typeof(t&time),getvardims(t&time))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
;filevardef(ncdf,"SKE",typeof(EKE),(/"time","lat","lon"/))
;filevardef(ncdf,"EKE",typeof(EKE),(/"time","lat","lon"/))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;define var attribute
;filevarattdef(ncdf,"time",t&time)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
;filevarattdef(ncdf,"EKE",EKE)
;filevarattdef(ncdf,"SKE",EKE)
filevarattdef(ncdf,"month_ave",EKE)
filevarattdef(ncdf,"season_ave",EKE)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
;ncdf->time = (/t&time/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
;ncdf->EKE  = (/EKE/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end

