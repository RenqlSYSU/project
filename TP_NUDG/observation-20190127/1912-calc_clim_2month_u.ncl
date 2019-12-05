;***************************************************************************
;                                            by Ql_Ren
;                                           2019/01/28
;******************************************************************************
begin
nc = 1
case = (/"ERA-Interim","NCEP1"/)
if(nc.eq.0) then 
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
else
path1   = "/home/ys17-19/data-observation/NCEP_NCAR/daily/"
end if
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd","hgt"/)
nv  = 0
nmonth = 4
month_s = 1
month_e = 8

years   = ispan(1979,2005,1)
years@calendar = "standard"

;area: east Asia
lats = -2
latn = 90
lonl =  0
lonr = 360
lev  = (/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100/) 

fileout  = path2 + "mdata/" + case(nc) + "_Clim_2month_" + var_name(nv) + ".nc"
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
if(nc.eq.0) then
f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
else
f    = addfile(path1+var_name(0)+"/"+var_name(0)+"."+years(0)+".nc", "r")
end if
vars = f->$var_name(0)$(0,0,{lats:latn},{lonl:lonr})
printVarSummary(vars)
print(vars(20:30,50:52))

time1 = cd_calendar(f->time,0)
select_time1 = ind(time1(:,1).ge.month_s.and.time1(:,1).le.month_e)
ntime = dimsizes(select_time1)

nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
var        = new((/nyear,ntime ,nlev,nlat,nlon/),float)
var_month  = new((/nyear,nmonth,nlev,nlat,nlon/),float)
month_ave  = new((/      nmonth,nlev,nlat,nlon/),float)
month_vari = new((/      nmonth,nlev,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups

;==========================================================
;read the data
;============================================================
if(ntime.gt.58) then 
    select_time2 = select_time1
    select_time2(58::) = select_time1(58::) + 1
end if

do nt = 0,nyear-1,1
    if(nc.eq.0) then
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    else
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+"."+years(nt)+".nc","r")
    end if
    if(isleapyear(years(nt))) then ;leap year,366 day
    var(nt,:,:,:,:) = f->$var_name(nv)$(select_time2,{lev},{lats:latn},{lonl:lonr})
    else
    var(nt,:,:,:,:) = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    end if
end do

do nm = 0,nmonth-1,1
    month_s1 = 2*nm + 1
    month_e1 = 2*nm + 2
    select_time := ind(time1(:,1).ge.month_s1.and.time1(:,1).le.month_e1)
    var_month(:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
end do

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
printVarSummary(month_ave)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "ERA-Interim uwnd from 1979-2005"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","lev","lat","lon"/);
dimSizes = (/nmonth ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"lev" ,typeof(var&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon" ,typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",var&level)
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"month_vari",var)

ncdf->lev  = (/lev/)
ncdf->lat  = (/var&lat/)
ncdf->lon  = (/var&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end

