;***************************************************************************
;filename:      calc_clim_month_ave_wind.ncl
;output_file:   Clim_month_ave_wind.nc 
;1.read the 27 years variable(TLLL) of one case 
;2.calculate the month average to get an array of (nmonth,nlat,nlon),
;3.calculate the season average (nseason,nlat,nlon).
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
    case    = (/"F2000","JET_NUDG","INT_NUDG"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_JET_NUDG/F2000_NUDG.cam.h1.YEAR.00",\
                      "/model/F2000_JET_NUDG2/F2000_NUDG2.cam.h1.YEAR.00"/)
    years    = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/2022EAJS/data/" + case + "-Clim_month_ave_"
end if
var_name = (/"U"  ,"V"    ,"T"  ,"Z3" ,"OMEGA","dudt","div",\
             "Q1" ,"preci","TMQ","TS" ,"PS"  , \
             "PSL","PRECC","PRECL" /) ;,"OMEGA","RELHUM" 

;Q1 is vertically integrated 

;area: east Asia
lats = 0  ;-90 ;;
latn = 80 ; 90  ;
lonl = 50 ; 0   ;
lonr = 180; 360 ;
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,150,100/) ;18 levels ,hPa U
;lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200/) ;hPa
;lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200,150,100,50,20, 10,5/) ;17levels hPa
;lev  = (/1000,925,850,700,500,400, 300,200,100,50,10/) ;,150,100hPa
;lev  = (/850,500,200/) ;v,

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
ntime = 365
var        = new((/nyear,365,nlev,nlat,nlon/),float)
var_month  = new((/nyear,12 ,nlev,nlat,nlon/),float)
var_season = new((/nyear,4  ,nlev,nlat,nlon/),float)
month_ave  = new((/      12 ,nlev,nlat,nlon/),float)
month_vari = new((/      12 ,nlev,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave = new((/       4 ,nlev,nlat,nlon/),float)
season_vari= new((/       4 ,nlev,nlat,nlon/),float)

do nv = 0,4,1; ,nvar-1 ,1 
do nf = 2,ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
if((nv.le.6)) then 
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
end if
if(nv.eq.6) then 
    u   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    v   = new((/nyear,ntime,nlev,nlat,nlon/),float)
    
    do nt = 0,nyear-1,1
    files := systemfunc("ls "+filein(nf)+"*.daily.U.nc" )
    print("handle with " + files(nt))
    f      = addfile(files(nt),"r")
    u(nt,:,:,:,:) = f->U(:,{lev},{lats:latn},{lonl:lonr})
    
    files := systemfunc("ls "+filein(nf)+"*.daily.V.nc" )
    print("handle with " + files(nt))
    f      = addfile(files(nt),"r")
    v(nt,:,:,:,:) = f->V(:,{lev},{lats:latn},{lonl:lonr})
    end do

    u = linmsg_n(u,-1,4)
    v = linmsg_n(v,-1,4)
    var = -uv2dv_cfd ( u, v, u&lat, u&lon, 2 )
    var@long_name = "divergence"
    var@units     = "1/s"
    delete([/u,v/])
end if
if(nv.le.4) then
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    do nt = 0,nyear-1,1
    print("handle with " + files(nt))
    f  = addfile(files(nt),"r")
;    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nt,:,:,:,:) = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    end do 
end if

time = cd_calendar(f->time,0)
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
end if

if(nv.ge.7) then
var        := new((/nyear,365,nlat,nlon/),float)
var_month  := new((/nyear,12 ,nlat,nlon/),float)
var_season := new((/nyear,4  ,nlat,nlon/),float)
month_ave  := new((/      12 ,nlat,nlon/),float)
month_vari := new((/      12 ,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups
season_ave := new((/       4 ,nlat,nlon/),float)
season_vari:= new((/       4 ,nlat,nlon/),float)
if(nv.eq.7) then 
var_name(nv) = "T"
files2 := systemfunc("ls "+filein(nf)+"*.daily.PS.nc" )
q1    = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
do nt = 0,nyear-1,1
    print("handle with " + files(nt))
    f  = addfile(files(nt),"r")
    q1(nt,:,:,:,:)  = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    f  = addfile(files2(nt),"r")
    ps(nt,:,:,:)   = f->PS(:,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
print(q1(10,10,2,10:20,30))
g  = 9.8 ;m/(s*s)
var = dim_sum_n(q1(:,:,5:10,:,:)*dp(:,:,5:10,:,:),2)/dim_sum_n(dp(:,:,5:10,:,:),2)
;var = dim_sum_n(q1*dp,2)/g
delete([/q1,ps,dp,g/])
;var@long_name = "Integrated Apparent Heat Source calc from temperature equation"
;var@units     = "W/kg"
var@long_name = "500-200hPa T"
var@units     = "K"
var_name(nv) = "uppT"
end if

if(nv.eq.8) then
    var = 0.0
    files1 := systemfunc("ls "+filein(nf)+"*.daily.PRECC.nc" )
    files2 := systemfunc("ls "+filein(nf)+"*.daily.PRECL.nc" )
    do nt = 0,nyear-1,1
    print("handle with " + files1(nt))
    f  = addfile(files1(nt),"r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->PRECC(:,{lats:latn},{lonl:lonr})
    f  = addfile(files2(nt),"r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->PRECL(:,{lats:latn},{lonl:lonr})
    end do 
    var  = var*1000*3600*24 ;convert m/s to mm/day
end if

if(nv.gt.8) then 
    files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
    do nt = 0,nyear-1,1
    print("handle with " + files(nt))
    f  = addfile(files(nt),"r")
    var(nt,:,:,:) = tofloat(f->$var_name(nv)$(:,{lats:latn},{lonl:lonr}))
    end do
    if (nv.eq.11.or.nv.eq.12) then
    var  = var*1000*3600*24 ;convert m/s to mm/day
    end if
end if

time = cd_calendar(f->time,0)
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
end if

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
season_ave  = dim_avg_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
season_vari = dim_variance_n_Wrap(var_season,0) ;get (nseason,nlat,nlon)
printVarSummary(month_ave)
printVarSummary(season_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+var_name(nv)+".nc")
ncdf = addfile(fileout(nf)+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the monthly and seasonal average"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"month","season","lev","lat","lon"/)
dimSizes = (/  12   ,  4     ,nlev ,nlat ,nlon/)
dimUnlim = (/False  ,False   ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
if((nv.le.6)) then 
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave",typeof(season_ave),(/"season","lev","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lev","lat","lon"/))
else
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave",typeof(season_ave),(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))
end if

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season  = (/1212,345,678,91011/)
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end do
end do
end

