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
var_name = (/"uwnd","vwnd","omega","air","hgt","ps"/)

years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats = -90 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
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

select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
ntime = dimsizes(select_time)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1

do nv = 5,5,1
year_option = True;False ;
fileout  = path2 + "ERA-Interim_DJF_" + var_name(nv) + "-year.nc"
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
    end do
else
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

;--------------------------------------
;DJF
;---------------------------------------------------
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
printVarSummary(var)
var1 = dim_avg_n_Wrap(var,1)
printVarSummary(var1)

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
end do
end

