;***************************************************************************
;filename:      calc_qtran_int_JJA.ncl
;use the post-processing data to calculate the integrated climatology qu, steady qu and 
;transient eddy qu, according to the Jiang Xinwen's paper 
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_JJAave_qtran_int_rmvannual.nc"
var_name= (/"U","V","Q","PS"/)
case    = (/"CTRL","NG24h","NG6h"/)
ncase   = dimsizes(case)
years   = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/1000,925,850,700,600,500,400,300/)
;lev  = (/600,500,400,300/)        ; high lev, according to Jiang Xinwen's paper
;lev  = (/1000,925,850,700,600/)   ;low lev

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime = dimsizes(time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(vars&level)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
q     = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)

do nf = 0,ncase-1,1 
;==========================================================
;read the data
;============================================================
do nt = 0,nyear-1,1
nv = 0
f      = addfile(filein(nf) + year(nt) + ".daily" + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)   = f->U(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
f      = addfile(filein(nf) + year(nt) + ".daily" + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->V(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
f      = addfile(filein(nf) + year(nt) + ".daily" + var_name(nv) + ".nc","r")
q(nt,:,:,:,:)   = f->Q(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
f      = addfile(filein(nf) + year(nt) + ".daily" + var_name(nv) + ".nc","r")
ps(nt,:,:,:)  = f->PS(select_time,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end do

printVarSummary(ps)
printVarSummary(u)
printVarSummary(q)

;calculates the pressure differences of a hybrid coordinate system.
dp  = dpres_plevel_Wrap(lev*100,ps,30000,0)  ;unit is Pa, get dp (year,time,lev,lat,lon)
printVarSummary(dp)

qu := dim_sum_n(q*u*dp,1)/g ;get integrated qu(nyear,time,lat,lon)
qv := dim_sum_n(q*v*dp,1)/g ;get integrated qu(nyear,time,lat,lon)
qu_vari = dim_variance_n(qu,0) 
qv_vari = dim_variance_n(qv,0) 
qu := dim_avg_n(qu,0)  ;get integrated qu(lat,lon)
qv := dim_avg_n(qv,0)  ;get integrated qu(lat,lon)

q_vari = dim_variance_n(dim_sum_n(q*dp,1),0)
u_vari = dim_variance_n(dim_sum_n(u*dp,1),0)
v_vari = dim_variance_n(dim_sum_n(v*dp,1),0)
q_ave = dim_avg_n(dim_sum_n(q*dp,1),0)
u_ave = dim_avg_n(dim_sum_n(u*dp,1),0)
v_ave = dim_avg_n(dim_sum_n(v*dp,1),0)
ste_qu = dim_sum_n(dim_avg_n(q,0)*dim_avg_n(u,0)*dim_avg_n(dp,0),0)/g ;get steady integrated qu(lat,lon)
ste_qv = dim_sum_n(dim_avg_n(q,0)*dim_avg_n(v,0)*dim_avg_n(dp,0),0)/g ;get steady integrated qv(lat,lon)

edd_qu := dim_sum_n(dim_rmvmean_n(q,0)*dim_rmvmean_n(u,0)*dp,1)/g ;get transient eddy integrated qu(lat,lon)
edd_qv := dim_sum_n(dim_rmvmean_n(q,0)*dim_rmvmean_n(v,0)*dp,1)/g
edd_qu_vari = dim_variance_n(edd_qu,0) ;get transient eddy integrated qu(lat,lon)
edd_qv_vari = dim_variance_n(edd_qv,0)
edd_qu := dim_avg_n(edd_qu,0) ;get transient eddy integrated qu(lat,lon)
edd_qv := dim_avg_n(edd_qv,0)

qu@long_name = "Integrated zonal moisture transport"
qv@long_name = "Integrated meridional moisture transport"
ste_qu@long_name = "steady Integrated meridional moisture transport"
edd_qu@long_name = "transient eddy Integrated meridional moisture transport"
qu@units     = "kg/(s*m)" 
qv@units     = "kg/(s*m)" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "the integration is based on 1000,925,850,700,500 hPa "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"lat" ,"lon"/)
dimSizes = (/nlat  ,nlon /)
dimUnlim = (/False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
filevardef(ncdf,"qu",typeof(qu),(/"lat","lon"/))
filevardef(ncdf,"qv",typeof(qv),(/"lat","lon"/))
filevardef(ncdf,"qu_vari",typeof(qu_vari),(/"lat","lon"/))
filevardef(ncdf,"qv_vari",typeof(qv_vari),(/"lat","lon"/))
filevardef(ncdf,"ste_qu",typeof(ste_qu),(/"lat","lon"/))
filevardef(ncdf,"ste_qv",typeof(ste_qv),(/"lat","lon"/))
filevardef(ncdf,"edd_qu",typeof(edd_qu),(/"lat","lon"/))
filevardef(ncdf,"edd_qv",typeof(edd_qv),(/"lat","lon"/))
filevardef(ncdf,"edd_qu_vari",typeof(edd_qu_vari),(/"lat","lon"/))
filevardef(ncdf,"edd_qv_vari",typeof(edd_qv_vari),(/"lat","lon"/))
filevardef(ncdf,"q_vari",typeof(q_vari),(/"lat","lon"/))
filevardef(ncdf,"u_vari",typeof(u_vari),(/"lat","lon"/))
filevardef(ncdf,"v_vari",typeof(v_vari),(/"lat","lon"/))
filevardef(ncdf,"q_ave",typeof(q_ave),(/"lat","lon"/))
filevardef(ncdf,"u_ave",typeof(u_ave),(/"lat","lon"/))
filevardef(ncdf,"v_ave",typeof(v_ave),(/"lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"qu",qu)
filevarattdef(ncdf,"qv",qv)
filevarattdef(ncdf,"edd_qu",edd_qu)
filevarattdef(ncdf,"edd_qv",edd_qv)

ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->qu   = (/qu/)
ncdf->qv   = (/qv/)
ncdf->ste_qu = (/ste_qu/)
ncdf->ste_qv = (/ste_qv/)
ncdf->edd_qu = (/edd_qu/)
ncdf->edd_qv = (/edd_qv/)
ncdf->qu_vari   = (/qu_vari/)
ncdf->qv_vari   = (/qv_vari/)
ncdf->edd_qu_vari = (/edd_qu_vari/)
ncdf->edd_qv_vari = (/edd_qv_vari/)
ncdf->q_vari = (/q_vari/)
ncdf->u_vari = (/u_vari/)
ncdf->v_vari = (/v_vari/)
ncdf->q_ave = (/q_ave/)
ncdf->u_ave = (/u_ave/)
ncdf->v_ave = (/v_ave/)
end do
end

