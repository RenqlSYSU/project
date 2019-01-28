;***************************************************************************
;                                            by Ql_Ren
;                                           2019/01/28
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd"/)
nv  = 0

years    = ispan(1979,2005,1)
years@calendar = "standard"

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/850,500,200/)
nl = 2

fileout  = path2 + "mdata/ERA-Interim_Clim_month_ave_" + lev(nl) + var_name(nv) + ".nc"
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,0,{lats:latn},{lonl:lonr})
ntime = 365
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
var   = new((/nyear,ntime,nlat,nlon/),float)
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

do nt = 0,nyear-1,1
if(isleapyear(years(nt))) then ;leap year,366 day
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    var(nt,:,:,:) = f->$var_name(nv)$(index,{lev(nl)},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    var(nt,:,:,:) = f->$var_name(nv)$(:,{lev(nl)},{lats:latn},{lonl:lonr})
end if
end do

f_model = addfile(path3,"r")
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
fileatt@discription   = "ERA-Interim uwnd from 1979-2005,850,500,200hPa"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","season","lat","lon"/);
dimSizes = (/  12   ,  4     ,nlat ,nlon/) ;
dimUnlim = (/False  ,False   ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"month" ,"integer",(/"month" /))
filevardef(ncdf,"season","integer",(/"season"/))
filevardef(ncdf,"lat" ,typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon" ,typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave)  ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_ave" ,typeof(season_ave) ,(/"season","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari) ,(/"month" ,"lat","lon"/))
filevardef(ncdf,"season_vari",typeof(season_vari),(/"season","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"season_ave",var)
filevarattdef(ncdf,"month_vari",var)
filevarattdef(ncdf,"season_vari",var)

ncdf->month  = ispan(1,12,1)
ncdf->season = (/1212,345,678,91011/)
ncdf->lat  = (/var&lat/)
ncdf->lon  = (/var&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->season_ave  = (/season_ave/)
ncdf->month_vari  = (/month_vari/)
ncdf->season_vari = (/season_vari/)
end

