;***************************************************************************
;                                            by Ql_Ren
;                                           2019/01/28
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","hgt","intQ1"/)
nv  = 3

years   = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats = 0  
latn = 90
lonl =  0
lonr = 360
;lev  = (/1000, 925, 850, 500, 200, 100/) 
lev  = (/1000,925,850,700,600,\
         500 ,400,350,300,250,\
         200 ,150,100/) ;hPa ;intQ1

do nv = 3,3,1
fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_" + var_name(nv) + ".nc"
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
;f    = addfile(path1+var_name(0)+"/"+var_name(0)+"."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
printVarSummary(vars)
print(vars(0,20:30,50:52))

ntime = 365
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)

;==========================================================
;read the data
;============================================================
index  = new((/ntime/),integer)
index1 = new((/ntime/),integer)
index(0:58)   = ispan(0,58,1)
index(59:364) = ispan(60,365,1)

f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)

if(nv.le.2) then 
var   = new((/nyear,ntime,nlev,nlat,nlon/),float)
var_month  = new((/nyear,12 ,nlev,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlev,nlat,nlon/),float)
month_ave  = new((/      12 ,nlev,nlat,nlon/),float)
month_vari = new((/      12 ,nlev,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlev,nlat,nlon/),float)
season_vari= new((/       4 ,nlev,nlat,nlon/),float)
do nt = 0,nyear-1,1
if(isleapyear(years(nt))) then ;leap year,366 day
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    ;f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    var(nt,:,:,:,:) = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    ;f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    var(nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
end if
end do

do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

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
end if

if(nv.eq.3) then  ;intQ1
var_month  := new((/nyear,12 ,nlat,nlon/),float)
var_season := new((/nyear,4  ,nlat,nlon/),float)
month_ave  := new((/      12 ,nlat,nlon/),float)
month_vari := new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave := new((/       4 ,nlat,nlon/),float)
season_vari:= new((/       4 ,nlat,nlon/),float)
u   = new((/nyear,ntime,nlev,nlat,nlon/),float)
v   = new((/nyear,ntime,nlev,nlat,nlon/),float)
t   = new((/nyear,ntime,nlev,nlat,nlon/),float)
w   = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp  = new((/nyear,ntime,nlev,nlat,nlon/),float)
f_sfcp  = addfile("/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc","r")
start = 0
do nt = 0,nyear-1,1
if(isleapyear(years(nt))) then ;leap year,366 day
    f      = addfile(path1+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:) = f->uwnd(index,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:) = f->vwnd(index,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"omega/omega.interim."+years(nt)+".nc","r")
    w(nt,:,:,:,:) = f->omega(index,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"air/air.interim."+years(nt)+".nc","r")
    t(nt,:,:,:,:) = f->air(index,{lev},{lats:latn},{lonl:lonr})
    space = 366
    index1(0:58)   = ispan(start,(start+58),1)
    index1(59:364) = ispan((start+60),(start+space-1),1)
    ps   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    f      = addfile(path1+"uwnd/uwnd.interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:) = f->uwnd(:,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"vwnd/vwnd.interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:) = f->vwnd(:,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"omega/omega.interim."+years(nt)+".nc","r")
    w(nt,:,:,:,:) = f->omega(:,{lev},{lats:latn},{lonl:lonr})
    f      = addfile(path1+"air/air.interim."+years(nt)+".nc","r")
    t(nt,:,:,:,:) = f->air(:,{lev},{lats:latn},{lonl:lonr})
    space = 365
    ps   = f_sfcp->sp(start:(start+space-1):1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps,20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
delete(ps)
printVarSummary(u)
printVarSummary(dp)

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
iday  = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
lev4d  = conform(u, lev*100, ilev)
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
loncyc = True  ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 
sigm0 = (R/cp)*(t/lev4d) - center_finite_diff_n(t,lev*100,cyclic,opt,ilev)
delete(lev4d)
coslat = conform(u, cos(lat),ilat)
Q1    =   (u/a/coslat) * center_finite_diff_n(t,lon,loncyc,opt,ilon) \
        + (v/a       ) * center_finite_diff_n(t,lat,cyclic,opt,ilat) \
        - sigm0*w + center_finite_diff_n(t,24*3600,cyclic,opt,iday)
delete([/u,v,t,w,coslat,sigm0/])
var = dim_sum_n(Q1*cp*dp,2)/g
delete(dp)
print(var(1,5,30,30))
var@long_name = "Apparent Heat Source calc from temperature equation"
var@units     = "W/kg" 

do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

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
end if

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
fileatt@discription   = "ERA-Interim uwnd from 1979-2005,850,500,200hPa"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","season","lev","lat","lon"/);
dimSizes = (/  12   ,  4     ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(nv.le.2) then 
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lev","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lev","lat","lon"/))
end if
if(nv.eq.3) then 
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do 
end

