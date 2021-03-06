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
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
nmonth = 2
month_s = 5
month_e = 8

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_2month_ave_advT.nc"
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_2month_ave_advT.nc"
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR"/) ;
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/energy-20180417/mdata/Ano2-Clim_2month_ave_qtran_int.nc"
end if

var_name= (/"U","V","Q","PS"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)

;area: east Asia
lats =  0  ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
;lev  = (/850,500,400,200/) 
lev  = (/1000,925,850,700,600,500,400,300,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+"01.daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)

;==========================================================
;read the data
;============================================================
qu_ave  = new((/2,nmonth,nlat,nlon/),double)
qv_ave  = new((/2,nmonth,nlat,nlon/),double)
div_ave = new((/2,nmonth,nlat,nlon/),double)
do nm = 0, nmonth-1, 1
mons = 2*nm + month_s
mone = 2*nm + month_s+1
print("month from " + mons + " to " + mone) 
select_time := ind(time(:,1).ge.mons.and.time(:,1).le.mone)
ntime = dimsizes(select_time)
u  := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
v  := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
q  := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
dp := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
ps := new((/ncase,nyear,ntime,nlat,nlon/),float)

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
do nt = 0,nyear-1,1
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
ps(nf,nt,:,:,:)  = f->$var_name(nv)$(select_time,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
dp(nf,nt,:,:,:,:)  = dpres_plevel_Wrap(lev*100,ps(nf,nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
end do
delete(ps)
printVarSummary(u)
u_month = dim_avg_n_Wrap(u,(/1,2/)) ;get (nyear,nmonth,nlat,nlon)
v_month = dim_avg_n_Wrap(v,(/1,2/)) ;get (nyear,nmonth,nlat,nlon)
q_month = dim_avg_n_Wrap(q,(/1,2/)) ;get (nyear,nmonth,nlat,nlon)
dp_month = dim_avg_n_Wrap(dp,(/1,2/)) ;get (nyear,nmonth,nlat,nlon)
delete([/v,q,dp/])

div_ave(0,nm,:,:)= dim_sum_n(-uv2dv_cfd((q_month(1,:,:,:)-q_month(0,:,:,:))*u_month(0,:,:,:), \ ;dim_avg_n(u_month,0)
                                        (q_month(1,:,:,:)-q_month(0,:,:,:))*v_month(0,:,:,:), \ ;dim_avg_n(v_month,0)
                                        u&lat, u&lon, 2)*dp_month(0,:,:,:),0)/g
div_ave(1,nm,:,:)= dim_sum_n(-uv2dv_cfd((u_month(1,:,:,:)-u_month(0,:,:,:))*q_month(0,:,:,:), \ ;dim_avg_n(q_month,0)
                                        (v_month(1,:,:,:)-v_month(0,:,:,:))*q_month(0,:,:,:), \ ;dim_avg_n(q_month,0)
                                        u&lat, u&lon, 2)*dp_month(0,:,:,:),0)/g
delete(u)
qu_ave(0,nm,:,:) = dim_sum_n((q_month(1,:,:,:)-q_month(0,:,:,:))*u_month(0,:,:,:)*dp_month(0,:,:,:),0)/g ;dim_avg_n(u_month,0)
qv_ave(0,nm,:,:) = dim_sum_n((q_month(1,:,:,:)-q_month(0,:,:,:))*v_month(0,:,:,:)*dp_month(0,:,:,:),0)/g ;dim_avg_n(v_month,0)
qu_ave(1,nm,:,:) = dim_sum_n((u_month(1,:,:,:)-u_month(0,:,:,:))*q_month(0,:,:,:)*dp_month(0,:,:,:),0)/g ;dim_avg_n(q_month,0)
qv_ave(1,nm,:,:) = dim_sum_n((v_month(1,:,:,:)-v_month(0,:,:,:))*q_month(0,:,:,:)*dp_month(0,:,:,:),0)/g ;dim_avg_n(q_month,0)
delete([/u_month,v_month,q_month,dp_month/])
end do
qu_ave@long_name = "steady moisture advection difference: q'*u, q*u'"
qu_ave@units     = "kg/m/s"
printVarSummary(qu_ave)
print(qu_ave(0,1,10,10:20))

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "rmvannual seasonal average"
fileattdef(ncdf,fileatt)  

;define all coordinate`
dimNames = (/"nvar","month","lat","lon"/);
dimSizes = (/   2  ,nmonth ,nlat ,nlon/) ;
dimUnlim = (/False ,False  ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"qu_ave"  ,typeof(qu_ave) ,(/"nvar","month" ,"lat","lon"/))
filevardef(ncdf,"qv_ave"  ,typeof(qv_ave) ,(/"nvar","month" ,"lat","lon"/))
filevardef(ncdf,"div_ave" ,typeof(div_ave),(/"nvar","month" ,"lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"qu_ave",qu_ave)
filevarattdef(ncdf,"qv_ave",qu_ave)
filevarattdef(ncdf,"div_ave",qu_ave)

ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->qu_ave   = (/qu_ave/)
ncdf->qv_ave   = (/qv_ave/)
ncdf->div_ave   = (/div_ave/)
end

