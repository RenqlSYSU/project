;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
;path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
path2   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","omega","air","hgt", \
             "ps"  ,"EGR","preci"/)

years    = ispan(1979,2016,1)
years@calendar = "standard"
DJF_option = False

;area: east Asia
lats =  5 ;-20
latn =  80 ;60
lonl =  5  ;20
lonr = 155 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
mdays = (/31,28,31,30,31,30,31,31,30,31,30,31/)

if(DJF_option) then
select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
ntime = dimsizes(select_time)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1
else
ntime = 365
select_time1  = new((/ntime/),integer)
select_time   = new((/ntime/),integer)
select_time1(0:58)   = ispan(0,58,1)
select_time1(59:364) = ispan(60,365,1)
select_time          = ispan(0,364,1)
end if

do nv = 6,6,1
if(DJF_option) then 
fileout  = path2 + "ERA-Interim_DJF_" + var_name(nv) + "-year.nc"
else
fileout  = path2 + "ERA-Interim_" + var_name(nv) + "-year.nc"
end if
;==========================================================
;read the data
;============================================================
if(nv.eq.5) then 
    var   := new((/nyear,ntime,nlat,nlon/),float)
    f = addfile("/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc","r")
    start = 0
    do nt = 0,nyear-1,1
    print("read: " + years(nt))
    if(isleapyear(years(nt))) then ;leap year,366 day
        space = 366
        index1 = select_time1 + start
        var(nt,:,:,:)   = f->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    else
        space = 365
        index1 = select_time + start
        var(nt,:,:,:)   = f->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
    end if
    start = start + space
    end do
end if

if(nv.lt.5) then 
    var   := new((/nyear,ntime,nlev,nlat,nlon/),float)
    do nt = 0,nyear-1,1
    print("read: " + years(nt))
    if(isleapyear(years(nt))) then ;leap year,366 day
        f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
        var(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    else   ;non-leap year, 365 day
        f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
        var(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    end if
    end do
end if

;------calc Eady growth rate--------------------------------
if(nv.eq.6) then 
    var   := new((/nyear,ntime,nlev,nlat,nlon/),float)
    air   := new((/ntime,nlev,nlat,nlon/),float)
    lat = conform(air,vars&lat,2) 
    
    lev_dim = 1
    opt = 0 ;opt=0, Return the Eady growth rate
;    opt=1, Return the Eady growth rate and the vertical gradient of the zonal wind (du/dz)
;    opt=2, Return the Eady growth rate and the vertical gradient of the zonal wind (du/dz) and the Brunt-Vaisala frequency

    opt    = 0     ;used by center_finite_diff_n, no meanging 
    cyclic = False ;used by center_finite_diff_n
    pi = atan(1.0)*4
    f0 = 2.0*(2.0*pi/24.0/3600.0)*sin(lat*pi/180.0)
    
    do nt = 0,nyear-1,1
    print("read: " + years(nt))
    if(isleapyear(years(nt))) then ;leap year,366 day
        f    = addfile(path1+"air/air.interim."+years(nt)+".nc","r")
        air  = f->air(select_time1,{lev},{lats:latn},{lonl:lonr})
        f    = addfile(path1+"hgt/hgt.interim."+years(nt)+".nc","r")
        hgt  = f->hgt(select_time1,{lev},{lats:latn},{lonl:lonr})
        f    = addfile(path1+"uwnd/uwnd.interim."+years(nt)+".nc","r")
        uwnd = f->uwnd(select_time1,{lev},{lats:latn},{lonl:lonr})
    else   ;non-leap year, 365 day
        f    = addfile(path1+"air/air.interim."+years(nt)+".nc","r")
        air  = f->air(select_time,{lev},{lats:latn},{lonl:lonr})
        f    = addfile(path1+"hgt/hgt.interim."+years(nt)+".nc","r")
        hgt  = f->hgt(select_time,{lev},{lats:latn},{lonl:lonr})
        f    = addfile(path1+"uwnd/uwnd.interim."+years(nt)+".nc","r")
        uwnd = f->uwnd(select_time,{lev},{lats:latn},{lonl:lonr})
    end if
    hgt = hgt/9.8 ;convert unit from m2/s2 to m
    
    ;static = static_stability(lev*100, air, lev_dim, 0)
    ;static = conform(air,dim_avg_n(static,0),(/1,2,3/))
    ;static = where(abs(static).lt.0.0000000001, 0.000001, static)
    ;var(nt,:,:,:,:)  = 0.3098*(f0/static)*center_finite_diff_n(uwnd, hgt, cyclic, opt, lev_dim)
    
    th  = pot_temp(lev*100, air, lev_dim, False)
    ;brunt = brunt_vaisala_atm(th, hgt, 0, lev_dim)
    ;brunt@_FillValue = -999
    ;brunt = where(brunt.le.0,brunt@_FillValue,brunt) 
    ;brunt = where(brunt.lt.0.00000000001,0.000001,brunt) 
    ;var(nt,:,:,:,:)  = 0.3098*(abs(f0)/brunt)*abs(center_finite_diff_n(uwnd, hgt, cyclic, opt, lev_dim))
    
    var(nt,:,:,:,:) = eady_growth_rate(th, uwnd, hgt, lat, opt, lev_dim)
    end do
    ;printVarSummary(brunt)
    ;printMinMax(brunt,True)
    printVarSummary(var)
    printMinMax(var,True)
    ;delete([/brunt,th,air,uwnd,hgt,lat/])
    ;delete([/static,air,uwnd,hgt,lat/])
    delete([/th,air,uwnd,hgt,lat/])
end if

;------------------------------------------------------
;calc month and season ave or DJF average
;---------------------------------------------------
if(DJF_option) then 
if(nv.eq.5) then 
    u2     := new((/nyear-1,ntime,nlat,nlon/),float)
    do ny = 0,nyear-2,1
        u2(ny,0:30,:,:) = var(ny,59:89,:,:)
        u2(ny,31:89,:,:) = var(ny+1,0:58,:,:)
    end do
else
    u2     := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
    do ny = 0,nyear-2,1
        u2(ny,0:30,:,:,:) = var(ny,59:89,:,:,:)
        u2(ny,31:89,:,:,:) = var(ny+1,0:58,:,:,:)
    end do
end if
var := u2
var1 = dim_avg_n_Wrap(var,1)

else

if(nv.eq.5) then 
var_month  := new((/nyear,12 ,nlat,nlon/),float)
var_season := new((/nyear,4  ,nlat,nlon/),float)
do nm = 0,11,1
    select_time2 := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:) = dim_avg_n_Wrap(var(:,select_time2,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do
;----------calculate the season average-----------------------
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
else
var_month  := new((/nyear,12 ,nlev,nlat,nlon/),float)
var_season := new((/nyear,4  ,nlev,nlat,nlon/),float)
do nm = 0,11,1
    select_time2 := ind(time(:,1).eq.nm+1)
    var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time2,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do
delete(var)
;----------calculate the season average-----------------------
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
end if

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "DJF mean" 
fileattdef(ncdf,fileatt)  

;define all coordinate
if(DJF_option) then
dimNames = (/"year" ,"lev","lat","lon"/);
dimSizes = (/nyear-1,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(nv.eq.5) then 
filevardef(ncdf,"season_ave",typeof(var1),(/"year","lat","lon"/))
else
filevardef(ncdf,"season_ave",typeof(var1),(/"year","lev","lat","lon"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"season_ave",var1)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->season_ave  = (/var1/)
else

dimNames = (/"year","month","season","lev","lat","lon"/)
dimSizes = (/nyear ,  12   ,  4     ,nlev ,nlat ,nlon/)
dimUnlim = (/False ,False  ,False   ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(nv.eq.5) then 
filevardef(ncdf,"month_ave" ,typeof(var_month) ,(/"year","month" ,"lat","lon"/))
filevardef(ncdf,"season_ave",typeof(var_season),(/"year","season","lat","lon"/))
else
filevardef(ncdf,"month_ave" ,typeof(var_month) ,(/"year","month" ,"lev","lat","lon"/))
filevardef(ncdf,"season_ave",typeof(var_season),(/"year","season","lev","lat","lon"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",vars)
filevarattdef(ncdf,"season_ave",vars)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->month_ave   = (/var_month/)
ncdf->season_ave  = (/var_season/)
end if
end do
end

