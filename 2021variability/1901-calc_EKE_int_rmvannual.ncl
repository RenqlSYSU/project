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
case_option = 3 ;0 means NUDG experiment; 1 means noTP experiment
option   = 1 ;0 means remove annual cycle, 1 means bandpass filter
calc_var = (/"intEKE","preci"/)
cnv = 1

;parameter for bw_bandpass_filter
ca = 80.0;90.0 ;day
cb = 8.0 ;2.5  ;10.0 ;

if(case_option.eq.3) then
    case    = (/"CTRL","NG15S","NG4550","NGEN","NG45N","NG1520"/) ;
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                      "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                      "/model/AMIP_NGEN/AMIP_CAM5_NGEN.cam.h1.YEAR.",\
                      "/model/AMIP_NG45N/AMIP_CAM5_NG45N.cam.h1.YEAR.",\
                      "/model/AMIP_NG1520/AMIP_CAM5_NG1520.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    if(option.eq.0) then
        fileout = path + "/project/2021variability/data/" + case + "-Clim_month_ave_"+calc_var(cnv)+".nc"
    else
        fileout = path + "/project/2021variability/data/" + case + "-Clim_month_ave_"+calc_var(cnv)+"-" + cb + "-" + ca + "day.nc"
    end if  
end if

ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -10   ;-90 ;
latn = 90  ; 90 ;
lonl = 0   ; 0  ;
lonr = 360 ;360 ;
lev  = (/1000,850,600,500,400, 300,200,100,50,10/) ;,150,100hPa
;lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100,50,20, 10,5/) 

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily.U.nc","r")
vars:= f->U(0,{lev(0)},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
ncase = dimsizes(filein)
var_month  = new((/nyear,12 ,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlat,nlon/),float)
month_ave  = new((/      12 ,nlat,nlon/),float)
month_vari = new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlat,nlon/),float)
season_vari= new((/       4 ,nlat,nlon/),float)

if(cnv.eq.0) then 
var_name= (/"U","V","PS"/)
u  = new((/nyear,ntime,nlev,nlat,nlon/),float)
v  = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps = new((/nyear,ntime,nlat,nlon/),float)
else
var_name= (/"PRECC","PRECL"/)
var = new((/nyear,ntime,nlat,nlon/),float)
end if

do nf = 0,3,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
if(cnv.eq.1) then 
do nt = 0,nyear-1,1
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
var(nt,:,:,:) = f->$var_name(nv)$(select_time,{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(select_time,{lats:latn},{lonl:lonr})
end do
var = linmsg_n(var,(/-1/),2)
var = var*3600*24*1000 ;mm/day
print(var(0,10,10:20,30))

fca = 1.0/ca
fcb = 1.0/cb
opt = False
var = bw_bandpass_filter(var,fca,fcb,opt,1)
var = var*var 
var@long_name = "variance"
var@units     = "mm2/day2" 
end if

if(cnv.eq.0) then 
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
filt_u = bw_bandpass_filter(u,fca,fcb,opt,0)
filt_v = bw_bandpass_filter(v,fca,fcb,opt,0)
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
cf = ezfftf_n(dim_avg_n(u,0),0)
cf(:,4:,:,:,:)  = 0
clm_u = ezfftb_n(cf,cf@xbar,1)
cf = ezfftf_n(dim_avg_n(v,0),0)
cf(:,4:,:,:,:)  = 0
clm_v = ezfftb_n(cf,cf@xbar,1)
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
var@long_name = "Integrated transient eddies kinetic energy";"Integrated stationary waves kinetic energy"
var@units     = "J" 
end if

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
fileatt@discription   = "" 
fileattdef(ncdf,fileatt)  

;define all coordinate`
dimNames = (/"month","season","lat","lon"/);
dimSizes = (/  12   ,  4     ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end

