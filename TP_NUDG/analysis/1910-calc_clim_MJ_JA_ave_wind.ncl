;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;read the 27 years variable(TLLL) of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the month average to get an array of (nmonth,nlat,nlon),
;then calculate the season average (nseason,nlat,nlon).
;
;Modification:
;20190712: add the case of noTP
;                                            by Ql_Ren
;                                           2018/01/18
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
var_name = (/"U","V","T","Z3","Q","dzdt","Ahgt","dTdy","dudt"/) ;,"OMEGA","RELHUM"
nv = 1
nmonth = 4
month_s = 1
month_e = 8

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_2month_ave_"
    years    = ispan(1979,2005,1)   ;the calculate time
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_"
    years   = ispan(21,30,1)   ;the calculate time
end if

if(case_option.eq.2) then
    case    = (/"F2000","NGall","NG58"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG5-8/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,2,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_"
end if
;area: east Asia
lats = 0 ;-90 ;
latn = 90  ;90  ;
lonl = 0  ;0   ;
lonr = 360 ;360 ;
lev  = (/1000,925,850,700,600,500,400,350,300,250,200,150,100/) ;hPa
;lev  = (/850,500,300,200,100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time1 = cd_calendar(f->time,0)
select_time1 = ind(time1(:,1).ge.month_s.and.time1(:,1).le.month_e)
ntime = dimsizes(select_time1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var        = new((/nyear,ntime ,nlev,nlat,nlon/),float)
var_month  = new((/nyear,nmonth,nlev,nlat,nlon/),float)
month_ave  = new((/      nmonth,nlev,nlat,nlon/),float)
month_vari = new((/      nmonth,nlev,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups

do nv = 2 ,4 ,1 
do nf = 0 , ncase-1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if(nv.eq.8) then 
    var_name2 = (/"U","V","OMEGA"/)
    u   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    v   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    w   = new((/nyear,ntime,nlev,nlat,nlon/),float)
end if

do nt = 0,nyear-1,1
if((nv.eq.5).or.(nv.eq.6)) then
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(3)+".nc","r")
    var(nt,:,:,:,:) = f->$var_name(3)$(select_time1,{lev},{lats:latn},{lonl:lonr})
end if
if(nv.eq.7) then
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(2)+".nc","r")
    var(nt,:,:,:,:) = f->$var_name(2)$(select_time1,{lev},{lats:latn},{lonl:lonr})
end if
if(nv.eq.8) then
    nv2 = 0
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    u(nt,:,:,:,:) = f->$var_name2(nv2)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv2 = 1
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    v(nt,:,:,:,:) = f->$var_name2(nv2)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv2 = 2
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    w(nt,:,:,:,:) = f->$var_name2(nv2)$(select_time1,{lev},{lats:latn},{lonl:lonr})
end if
if(nv.lt.5) then
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f      = addfile(files(nt),"r")
    var(nt,:,:,:,:) = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
end if
end do

if(nv.eq.5) then ;calc dzdt
    opt    = 0     ;used by center_finite_diff_n, no meanging 
    cyclic = False ;used by center_finite_diff_n
    ;means NOT treated as cyclic in differenced coordinates
    ;the end values will use a one-sided difference scheme for the end points. 
    var = linmsg_n(var,-1,4)
    var = center_finite_diff_n(var,24*3600,cyclic,opt,1)
    var@long_name = "geopotential tendencyi tendency"
    var@units     = "gpm/s" 
end if

if(nv.eq.8) then 
    iyear = 0
    itime = 1
    ilev  = 2
    ilat  = 3
    ilon  = 4
    pi = atan(1.0)*4
    a  = 6378388 ;the radius of earth, m
    lat  = vars&lat
    lon  = vars&lon
    lat  = lat*pi/180.0
    lon  = lon*pi/180.0
    coslat = conform(u,cos(lat),ilat)
    opt    = 0     ;used by center_finite_diff_n, no meanging 
    loncyc = True ;used by center_finite_diff_n
    cyclic = False ;used by center_finite_diff_n
    term = (u/a/coslat)*center_finite_diff_n(u,lon,loncyc,opt,ilon) + \
           (v/a       )*center_finite_diff_n(u,lat,cyclic,opt,ilat) + \
           w*center_finite_diff_n(u,lev*100,cyclic,opt,ilev) + \
           center_finite_diff_n(u,24*3600,cyclic,opt,itime)
    var = tofloat(term)
    var@long_name = "dudt = dudt + u*dudx + v*dudy + w*dudz"
    var@units     = "m*s^-2"
    delete([/u,v,w,term,lat,lon,pi,coslat,a/])
end if

time = cd_calendar(f->time(select_time1),0)
do nm = 0, nmonth-1, 1
    month_s = 2*nm + 1
    month_e = 2*nm + 2
    select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
    var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

if(nv.eq.6) then ;calc dzdt
var_month = var_month - conform(var_month,dim_avg_n(var_month,4),(/0,1,2,3/))
var_month@long_name = "zonal asymmetric geopotential height"
var_month@units     = "gpm" 
end if

if(nv.eq.7) then ;calc dzdt
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
lat = var&lat
lat = lat*pi/180.0
coslat = conform(var_month,cos(lat),3)
var_month = tofloat(center_finite_diff_n(var_month*coslat,lat,False,0,3)/a/coslat)
var_month@long_name = "dTdy, meridional gradient of temperature"
var_month@units     = "K/m" 
end if

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
copy_VarMeta(var(0,0,:,:,:),var_month(0,0,:,:,:))
printVarSummary(month_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+var_name(nv)+".nc")
ncdf = addfile(fileout(nf)+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "ealy(56) and late(78) stage of the summer"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","lev","lat","lon"/)
dimSizes = (/nmonth ,nlev ,nlat ,nlon/)
dimUnlim = (/False  ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari",typeof(month_vari),(/"month" ,"lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var_month)
filevarattdef(ncdf,"month_vari",var_month)

ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end do
end do ;loop of variable
end

