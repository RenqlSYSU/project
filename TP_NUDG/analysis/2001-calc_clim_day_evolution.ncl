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

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_day_evolution_"
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
    case    = (/"F2000","NUDG"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_"
end if
var_name = (/"U","V","T","Z3","Q1","dudt"/) ;,"OMEGA","RELHUM"

;area: east Asia
lats = 0   ;-90 ;
latn = 90  ;90  ;
lonl = 0   ;0   ;
lonr = 360 ;360 ;
lev  = (/500,200/) ;v,

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time = f->time
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
ntime = 365
var        = new((/nyear,365,nlev,nlat,nlon/),float)
var_day    = new((/      365,nlev,nlat,nlon/),float)

do nv = 0,0,1; ,nvar-1 ,1 
do nf = 0 , ncase-1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if((nv.le.3).or.(nv.eq.5)) then 
if(nv.eq.5) then 
    var_name2 = (/"U","V","OMEGA"/)
    u   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    v   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    w   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    
    do nt = 0,nyear-1,1
    nv2 = 0
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    u(nt,:,:,:,:) = f->$var_name2(nv2)$(:,{lev},{lats:latn},{lonl:lonr})
    nv2 = 1
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    v(nt,:,:,:,:) = f->$var_name2(nv2)$(:,{lev},{lats:latn},{lonl:lonr})
    nv2 = 2
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name2(nv2)+".nc" )
    f      = addfile(files(nt),"r")
    w(nt,:,:,:,:) = f->$var_name2(nv2)$(:,{lev},{lats:latn},{lonl:lonr})
    end do

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
else
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    do nt = 0,nyear-1,1
    print("handle with " + years(nt))
    f  = addfile(files(nt),"r")
;    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nt,:,:,:,:) = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    end do 
end if
end if

if(nv.eq.4) then 
var        := new((/nyear,365,nlat,nlon/),float)
var_day    := new((/      365,nlat,nlon/),float)
files2 := systemfunc("ls "+filein(nf)+"*.daily.PS.nc" )
q1    = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)
do nt = 0,nyear-1,1
    print("handle with " + years(nt))
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    f  = addfile(files(nt),"r")
    q1(nt,:,:,:,:)  = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    f  = addfile(files2(nt),"r")
    ps(nt,:,:,:)   = f->PS(:,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
print(q1(10,10,2,10:20,30))
g  = 9.8 ;m/(s*s)
var = dim_sum_n(q1*dp,2)/g
delete([/q1,ps,dp,g/])
var@long_name = "Integrated Apparent Heat Source calc from temperature equation"
var@units     = "W/kg" 
end if

var_day = dim_avg_n_Wrap(var,0)
printVarSummary(var_day)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+var_name(nv)+".nc")
ncdf = addfile(fileout(nf)+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the daily average"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"time","lev","lat","lon"/)
dimSizes = (/ntime ,nlev ,nlat ,nlon/)
dimUnlim = (/False ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(time),getvardims(time))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
if((nv.le.3).or.(nv.eq.5)) then 
filevardef(ncdf,"var_day" ,typeof(var_day) ,(/"time" ,"lev","lat","lon"/))
else
filevardef(ncdf,"var_day" ,typeof(var_day) ,(/"time" ,"lat","lon"/))
end if

;Define the attribute of the variables
filevarattdef(ncdf,"time",time)
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"var_day",var_day)

ncdf->time = (/time/)
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->var_day = (/var_day/)
end do
end do
end

