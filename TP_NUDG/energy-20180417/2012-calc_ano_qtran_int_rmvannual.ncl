;***************************************************************************
;use the post-processing data to calculate the integrated climatology qu, steady qu and 
;transient eddy qu, and this transient eddy is calc based an annual cycle 
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/energy-20180417/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_qtran_int_rmvannual.nc"
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/energy-20180417/mdata/" + case + "-Clim_qtran_int_rmvannual.nc"
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR"/) ;,"NUDG" 
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00"/)
;                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/energy-20180417/mdata/Ano2-Clim_qtran_int_rmvannual.nc"
end if

var_name= (/"U","V","Q","PS"/)
;var_name= (/"U","V","T","PS"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)

;area: east Asia
lats =  10  ;-20
latn =  70 ;60
lonl =  20 ;20
lonr = 160 ;220
lev  = (/1000,925,850,700,600,500,400,300,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+"01.daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
u     = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
q     = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/ncase,nyear,ntime,nlat,nlon/),float)

;==========================================================
;read the data
;============================================================
do nf = 0,1,1 
do nt = 0,nyear-1,1
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nf,nt,:,:,:,:)   = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q(nf,nt,:,:,:,:)   = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
ps(nf,nt,:,:,:)  = f->$var_name(nv)$(select_time,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
dp(nf,nt,:,:,:,:)  = dpres_plevel_Wrap(lev*100,ps(nf,nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
end do
u = linmsg_n(u,(/-1/),5)
v = linmsg_n(v,(/-1/),5)
q = linmsg_n(q,(/-1/),5)

printVarSummary(ps)
printVarSummary(u)
printVarSummary(q)
printVarSummary(dp)

if(case_option.ge.1) then 
    cf = ezfftf_n(dim_avg_n(u,1),1)
    cf(:,:,4:,:,:,:)  = 0
    clm_u = ezfftb_n(cf,cf@xbar,2)
    cf = ezfftf_n(dim_avg_n(v,1),1)
    cf(:,:,4:,:,:,:)  = 0
    clm_v = ezfftb_n(cf,cf@xbar,2)
    cf = ezfftf_n(dim_avg_n(q,1),1)
    cf(:,:,4:,:,:,:)  = 0
    clm_q = ezfftb_n(cf,cf@xbar,2)
    cf = ezfftf_n(dim_avg_n(dp,1),1)
    cf(:,:,4:,:,:,:)  = 0
    clm_dp = ezfftb_n(cf,cf@xbar,2)
else
    clm_u = dim_avg_n(u,1)
    clm_v = dim_avg_n(v,1)
    clm_q = dim_avg_n(q,1)
    clm_dp = dim_avg_n(dp,1)
end if
delete([/v,q,dp/])
div1 = dim_sum_n(-uv2dv_cfd((clm_q(1,:,:,:,:)-clm_q(0,:,:,:,:))*clm_u(0,:,:,:,:), \ ;dim_avg_n(clm_u,0)
                            (clm_q(1,:,:,:,:)-clm_q(0,:,:,:,:))*clm_v(0,:,:,:,:), \ ;dim_avg_n(clm_v,0)
                            u&lat, u&lon, 2)*clm_dp(0,:,:,:,:),1)/g
div2 = dim_sum_n(-uv2dv_cfd((clm_u(1,:,:,:,:)-clm_u(0,:,:,:,:))*clm_q(0,:,:,:,:), \ ;dim_avg_n(clm_q,0)
                            (clm_v(1,:,:,:,:)-clm_v(0,:,:,:,:))*clm_q(0,:,:,:,:), \ ;dim_avg_n(clm_q,0)
                            u&lat, u&lon, 2)*clm_dp(0,:,:,:,:),1)/g
delete(u)
qu1 = dim_sum_n((clm_q(1,:,:,:,:)-clm_q(0,:,:,:,:))*clm_u(0,:,:,:,:)*clm_dp(0,:,:,:,:),1)/g ;dim_avg_n(clm_u,0)
qv1 = dim_sum_n((clm_q(1,:,:,:,:)-clm_q(0,:,:,:,:))*clm_v(0,:,:,:,:)*clm_dp(0,:,:,:,:),1)/g ;dim_avg_n(clm_v,0)
qu2 = dim_sum_n((clm_u(1,:,:,:,:)-clm_u(0,:,:,:,:))*clm_q(0,:,:,:,:)*clm_dp(0,:,:,:,:),1)/g ;dim_avg_n(clm_q,0)
qv2 = dim_sum_n((clm_v(1,:,:,:,:)-clm_v(0,:,:,:,:))*clm_q(0,:,:,:,:)*clm_dp(0,:,:,:,:),1)/g ;dim_avg_n(clm_q,0)
delete([/clm_u,clm_v,clm_q,clm_dp/])

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "the integration is based on 1000,925,850,700,600,500,400,300,200 hPa "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lat" ,"lon"/)
dimSizes = (/ntime ,nlat  ,nlon /)
dimUnlim = (/False ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"qu1",typeof(qu1),(/"time","lat","lon"/))
filevardef(ncdf,"qv1",typeof(qv1),(/"time","lat","lon"/))
filevardef(ncdf,"qv2",typeof(qv2),(/"time","lat","lon"/))
filevardef(ncdf,"qu2",typeof(qu2),(/"time","lat","lon"/))
filevardef(ncdf,"div1",typeof(div1),(/"time","lat","lon"/))
filevardef(ncdf,"div2",typeof(div2),(/"time","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"qu1",qu1)

ncdf->time = (/vars&time/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->qu1   = (/qu1/)
ncdf->qv1   = (/qv1/)
ncdf->qu2   = (/qu2/)
ncdf->qv2   = (/qv2/)
ncdf->div1  = (/div1/)
ncdf->div2  = (/div2/)
end

