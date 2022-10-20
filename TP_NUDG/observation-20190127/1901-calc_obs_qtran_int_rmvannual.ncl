;***************************************************************************
;use the post-processing data to calculate the integrated climatology qu, steady qu and 
;transient eddy qu, and this transient eddy is calc based an annual cycle 
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc"
var_name = (/"uwnd","vwnd","shum"/)
fileout  = path2 + "mdata/Obs_Clim_daily_qtran.nc"
g  = 9.8 ;m/(s*s)
years    = ispan(1979,2005,1)
years@calendar = "standard"

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/1000,925,850,700,600,500,400,300/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
    t   = f_model->U(:,0,lats,lonl)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,0,{lats:latn},{lonl:lonr})
ntime = 365
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
u     = new((/nyear,ntime,nlev,nlat,nlon/),float)
v     = new((/nyear,ntime,nlev,nlat,nlon/),float)
q     = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp    = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps    = new((/nyear,ntime,nlat,nlon/),float)

;==========================================================
;read the data
;============================================================
index  = new((/ntime/),integer)
index1 = new((/ntime/),integer)
index(0:58)   = ispan(0,58,1)
index(59:364) = ispan(60,365,1)

f_sfcp  = addfile(path4,"r")
start = 0
do nt = 0,nyear-1,1
print("readi: " + years(nt) + ", start = " + start)
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    q(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    space = 366
    index1(0:58)   = ispan(start,(start+58),1)
    index1(59:364) = ispan((start+60),(start+space-1),1)
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    nv = 2
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    q(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    space = 365
    ps(nt,:,:,:)   = f_sfcp->sp(start:(start+space-1):1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),30000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
printVarSummary(ps)
printVarSummary(u)
printVarSummary(q)
printVarSummary(dp)

qu := dim_sum_n(q*u*dp,2)/g ;get integrated qu(nyear,time,lat,lon)
qv := dim_sum_n(q*v*dp,2)/g ;get integrated qu(nyear,time,lat,lon)
qu := dim_avg_n(qu,0)  ;get integrated qu(lat,lon)
qv := dim_avg_n(qv,0)  ;get integrated qu(lat,lon)

ste_qu := dim_sum_n(dim_avg_n(q,0)*dim_avg_n(u,0)*dim_avg_n(dp,0),1)/g ;get steady integrated qu(ntime,lat,lon)
ste_qv := dim_sum_n(dim_avg_n(q,0)*dim_avg_n(v,0)*dim_avg_n(dp,0),1)/g ;get steady integrated qv(ntime,lat,lon)

edd_q = q - conform(q,dim_avg_n(q,0),(/1,2,3,4/))
edd_u = u - conform(u,dim_avg_n(u,0),(/1,2,3,4/))
edd_v = v - conform(v,dim_avg_n(v,0),(/1,2,3,4/))
edd_qu := dim_sum_n(edd_q*edd_u*dp,2)/g ;get transient eddy integrated qu(nyear,ntime,lat,lon)
edd_qv := dim_sum_n(edd_q*edd_v*dp,2)/g
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
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "the integration is based on 1000,925,850,700,600,500,400,300 hPa "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lat" ,"lon"/)
dimSizes = (/ntime ,nlat  ,nlon /)
dimUnlim = (/False ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(t&time),getvardims(t&time))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
filevardef(ncdf,"qu",typeof(qu),(/"time","lat","lon"/))
filevardef(ncdf,"qv",typeof(qv),(/"time","lat","lon"/))
filevardef(ncdf,"ste_qu",typeof(ste_qu),(/"time","lat","lon"/))
filevardef(ncdf,"ste_qv",typeof(ste_qv),(/"time","lat","lon"/))
filevardef(ncdf,"edd_qu",typeof(edd_qu),(/"time","lat","lon"/))
filevardef(ncdf,"edd_qv",typeof(edd_qv),(/"time","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",t&time)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"qu",qu)
filevarattdef(ncdf,"qv",qv)
filevarattdef(ncdf,"edd_qu",edd_qu)
filevarattdef(ncdf,"edd_qv",edd_qv)

ncdf->time = (/t&time/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->qu   = (/qu/)
ncdf->qv   = (/qv/)
ncdf->ste_qu = (/ste_qu/)
ncdf->ste_qv = (/ste_qv/)
ncdf->edd_qu = (/edd_qu/)
ncdf->edd_qv = (/edd_qv/)
end

