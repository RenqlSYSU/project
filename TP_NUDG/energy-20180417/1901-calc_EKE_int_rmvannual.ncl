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
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5  ;10.0 ;

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_ave_500EKE.nc"
        ;fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_ave_intuv.nc"
    else
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_month_ave_intEKE2-" + cb + "-" + ca + "day.nc"
    end if  
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    if(option.eq.0) then
        ;fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_uppEKE.nc"
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_intEKE.nc"
    else
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_lowEKE-" + cb + "-" + ca + "day.nc"
    end if  
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR","NUDG"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_intEKE.nc"
        ;fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_intuv.nc"
    else
        fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_intEKE-" + cb + "-" + ca + "day.nc"
    end if  
end if

var_name= (/"U","V","PS"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = 0   ;-90 ;
latn = 90  ; 90 ;
lonl = 0   ; 0  ;
lonr = 360 ;360 ;
lev  = (/1000,925,850,700,600,500,400,350,300,250,200/) ;,150,100hPa

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
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

do nf = 1,1,1 
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nt = 0,nyear-1,1
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f   = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
ps(nt,:,:,:)   = f->$var_name(nv)$(select_time,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
u = linmsg_n(u,(/-1/),1)
v = linmsg_n(v,(/-1/),1)
print(u(5,10,2,10:20,30))

if(option.eq.1) then
u := reshape(u,(/nyear*ntime,nlev,nlat,nlon/))
v := reshape(v,(/nyear*ntime,nlev,nlat,nlon/))
fca = 1.0/ca
fcb = 1.0/cb
opt = False
filt_u = bw_bandpass_filter(u,fca,fcb,opt,1)
filt_v = bw_bandpass_filter(v,fca,fcb,opt,1)
delete([/u,v/])
filt_u := reshape(filt_u,(/nyear,ntime,nlev,nlat,nlon/))
filt_v := reshape(filt_v,(/nyear,ntime,nlev,nlat,nlon/))

;fcf = fspan(0, 0.5, ntime/2) 
;filt_cf = ind(fcf.lt.fca.or.fcf.gt.fcb)
;cf = ezfftf_n(u,1)
;cf(:,:,filt_cf,:,:,:) = 0.0  
;filt_u = ezfftb_n(cf,0,2)
;cf = ezfftf_n(v,1)
;cf(:,:,filt_cf,:,:,:) = 0.0  
;filt_v = ezfftb_n(cf,0,2)

;var   = dim_sum_n((filt_v*filt_v + filt_u*filt_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
;var   = dim_sum_n((filt_v(:,:,5:12,:,:)*filt_v(:,:,5:12,:,:) + filt_u(:,:,5:12,:,:)*filt_u(:,:,5:12,:,:))*dp(:,:,5:12,:,:),2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
var   = dim_sum_n((filt_v(:,:,0:5,:,:)*filt_v(:,:,0:5,:,:) + filt_u(:,:,0:5,:,:)*filt_u(:,:,0:5,:,:))*dp(:,:,0:5,:,:),2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
delete(filt_u)
delete(filt_v)
else
if(case_option.ge.1) then 
    cf = ezfftf_n(dim_avg_n(u,0),0)
    cf(:,4:,:,:,:)  = 0
    clm_u = ezfftb_n(cf,cf@xbar,1)
    cf = ezfftf_n(dim_avg_n(v,0),0)
    cf(:,4:,:,:,:)  = 0
    clm_v = ezfftb_n(cf,cf@xbar,1)
else
    clm_u = dim_avg_n(u,0)
    clm_v = dim_avg_n(v,0)
end if
edd_u = u - conform(u,clm_u,(/1,2,3,4/))
edd_v = v - conform(v,clm_v,(/1,2,3,4/))   
;ste_u = conform(u,dim_avg_n(u,0),(/1,2,3,4/))
;ste_v = conform(v,dim_avg_n(v,0),(/1,2,3,4/))
;var   = dim_sum_n((edd_u*edd_v)*dp,2)/g/2.0   ;ps-200hPa integral get (nyear,ntime,nlat,nlon)
var   = dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,2)/g/2.0   ;ps-200hPa integral get (nyear,ntime,nlat,nlon)
;var   = dim_sum_n( (edd_v(:,:,0:5,:,:)*edd_v(:,:,0:5,:,:) + edd_u(:,:,0:5,:,:)*edd_u(:,:,0:5,:,:))*dp(:,:,0:5,:,:), 2 )/g/2.0   ;surface-500hPa integralget (nyear,ntime,nlat,nlon)
;var   = dim_sum_n( (edd_v(:,:,5:10,:,:)*edd_v(:,:,5:10,:,:) + edd_u(:,:,5:10,:,:)*edd_u(:,:,5:10,:,:))*dp(:,:,5:10,:,:), 2 )/g/2.0   ;500-200hPa integralget (nyear,ntime,nlat,nlon)
;var   = dim_sum_n((ste_v*ste_v + ste_u*ste_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
;var   = dim_sum_n((ste_v(:,:,5:8,:,:)*ste_v(:,:,5:8,:,:) + ste_u(:,:,5:8,:,:)*ste_u(:,:,5:8,:,:))*dp(:,:,5:8,:,:),2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
delete([/edd_u,edd_v/])
end if

EKE   = dim_avg_n(var,0) 
;EKE_vari = dim_variance_n(var,0) 

;EKE@long_name = "Integrated eddy momentum flux" 
;EKE@units     = "m2/s2" 
EKE@long_name = "Integrated transient eddies kinetic energy";"Integrated stationary waves kinetic energy"
EKE@units     = "J" 

time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time1 := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time1,:,:),1) ;get (nyear,nmonth,nlat,nlon)
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
print(month_ave(6,10,10:20))

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "rmvannual cycle, the integration is based on 1000,925,850,700,600,500hPa 500,400,300,200"
fileattdef(ncdf,fileatt)  

;define all coordinate`
dimNames = (/"month","season","time","lat","lon"/);
dimSizes = (/  12   ,  4     ,ntime ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"time",typeof(u&time),getvardims(u&time))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
;filevardef(ncdf,"uv",typeof(EKE),(/"time","lat","lon"/))
;filevardef(ncdf,"uv_vari",typeof(EKE_vari),(/"time","lat","lon"/))
;filevardef(ncdf,"SKE",typeof(EKE),(/"time","lat","lon"/))
;filevardef(ncdf,"SKE_vari",typeof(EKE_vari),(/"time","lat","lon"/))
;filevardef(ncdf,"EKE",typeof(EKE),(/"time","lat","lon"/))
;filevardef(ncdf,"EKE_vari",typeof(EKE_vari),(/"time","lat","lon"/))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
;filevarattdef(ncdf,"uv",EKE)
;filevarattdef(ncdf,"EKE",EKE)
;filevarattdef(ncdf,"SKE",EKE)
filevarattdef(ncdf,"month_ave",EKE)
filevarattdef(ncdf,"season_ave",EKE)
filevarattdef(ncdf,"month_vari",EKE)
filevarattdef(ncdf,"season_vari",EKE)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->time = (/u&time/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
;ncdf->uv  = (/EKE/)
;ncdf->uv_vari  = (/EKE_vari/)
;ncdf->SKE  = (/EKE/)
;ncdf->SKE_vari  = (/EKE_vari/)
;ncdf->EKE  = (/EKE/)
;ncdf->EKE_vari  = (/EKE_vari/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end

