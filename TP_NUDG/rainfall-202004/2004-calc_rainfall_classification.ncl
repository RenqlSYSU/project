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
nmonth = 4
month_s = 1
month_e = 8

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/rainfall-202004/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-rainfall_frequency.nc"
    years    = ispan(1979,2005,1)   ;the calculate time
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_preci.nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

var_name = (/"PRECC","PRECL"/) ;,"OMEGA","RELHUM"

;area: east Asia
lats =  0  ;-90 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
time1 = cd_calendar(f->time,0)
select_time1 = ind(time1(:,1).ge.month_s.and.time1(:,1).le.month_e)
ntime = dimsizes(select_time1)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
ntype = 5
var        = new((/nyear,ntime ,nlat,nlon/),float)
var_month  = new((/nyear,ntype,nmonth,nlat,nlon/),float)
month_ave  = new((/      ntype,nmonth,nlat,nlon/),float)
month_vari = new((/      ntype,nmonth,nlat,nlon/),float) ;be used to exam the significance of the mean difference between the two groups

do nf = 0 , ncase-1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
var = 0
do nv = 0,nvar-1,1 
do nt = 0,nyear-1,1
    f   = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(select_time1,{lats:latn},{lonl:lonr})
end do
end do
var  = var*1000*3600*24 ;convert m/s to mm/day

time = cd_calendar(f->time(select_time1),0)
do nm = 0, nmonth-1, 1
mons = 2*nm + 1
mone = 2*nm + 2
select_time := ind(time(:,1).ge.mons.and.time(:,1).le.mone)
nday = dimsizes(select_time)
print("===================================================================")
print("month from " + mons + " to " + mone + ", ntime = " + nday)

term := where(var(:,select_time,:,:).lt.0.1,1.0,0.0) ;no rain
var_month(:,0,nm,:,:) = dim_sum_n_Wrap(term,1)/tofloat(nday) 
print("no rain")
print(max(var_month(:,0,nm,:,:)))
print(min(var_month(:,0,nm,:,:)))
term := where(var(:,select_time,:,:).lt.10.0.and.var(:,select_time,:,:).ge.0.1,1.0,0.0) ;small rain
var_month(:,1,nm,:,:) = dim_sum_n_Wrap(term,1)/tofloat(nday) 
print("small rain")
print(max(var_month(:,1,nm,:,:)))
print(min(var_month(:,1,nm,:,:)))
term := where(var(:,select_time,:,:).lt.25.0.and.var(:,select_time,:,:).ge.10.0,1.0,0.0) ;medium rain
var_month(:,2,nm,:,:) = dim_sum_n_Wrap(term,1)/tofloat(nday) 
print("medium rain")
print(max(var_month(:,2,nm,:,:)))
print(min(var_month(:,2,nm,:,:)))
term := where(var(:,select_time,:,:).lt.50.0.and.var(:,select_time,:,:).ge.25.0,1.0,0.0) ;large rain
var_month(:,3,nm,:,:) = dim_sum_n_Wrap(term,1)/tofloat(nday) 
print("large rain")
print(max(var_month(:,3,nm,:,:)))
print(min(var_month(:,3,nm,:,:)))
term := where(var(:,select_time,:,:).ge.50.0,1.0,0.0) ;storm
var_month(:,4,nm,:,:) = dim_sum_n_Wrap(term,1)/tofloat(nday) 
print("storm")
print(max(var_month(:,4,nm,:,:)))
print(min(var_month(:,4,nm,:,:)))
end do

month_ave  = dim_avg_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,0) ;get (nmonth,nlat,nlon)
printVarSummary(month_ave)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "winter(1,2), spring(3,4), ealy(56) and late(78) stage of the summer" +\
                        "norain, smallrain, mediumrain, largerain, storm"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"type","month","lat","lon"/)
dimSizes = (/ntype ,nmonth ,nlat ,nlon/)
dimUnlim = (/False ,False  ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave" ,typeof(month_ave) ,(/"type","month","lat","lon"/))
filevardef(ncdf,"month_vari",typeof(month_vari),(/"type","month","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"month_vari",var)

ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end do
;end do ;loop of variable
end

