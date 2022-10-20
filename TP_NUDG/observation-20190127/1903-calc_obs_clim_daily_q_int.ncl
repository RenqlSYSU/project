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
var_name = (/"shum"/)
fileout  = path2 + "mdata/ERA-Interim_Clim_daily_q_int.nc"
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
    t   = f_model->U(:,0,{lats},{lonl})

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,0,{lats:latn},{lonl:lonr})
ntime = 365
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
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
print("read: " + years(nt) + ", start = " + start)
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    q(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats:latn},{lonl:lonr})
    space = 366
    index1(0:58)   = ispan(start,(start+58),1)
    index1(59:364) = ispan((start+60),(start+space-1),1)
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    q(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    space = 365
    ps(nt,:,:,:)   = f_sfcp->sp(start:(start+space-1):1,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),30000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
printVarSummary(ps)
printVarSummary(q)
printVarSummary(dp)
q_int = dim_avg_n(dim_sum_n(q*dp,2)/g,0) ;get integrated qu(nyear,time,lat,lon)
q_int@long_name = "vertically integrated water vapor"
q_int@units     = "kg/m~S~2~N~" 

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
filevardef(ncdf,"lat" ,typeof(q&lat),getvardims(q&lat))
filevardef(ncdf,"lon" ,typeof(q&lon),getvardims(q&lon))
filevardef(ncdf,"q",typeof(q_int),(/"time","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",t&time)
filevarattdef(ncdf,"lat",q&lat)
filevarattdef(ncdf,"lon",q&lon)
filevarattdef(ncdf,"q",q_int)

ncdf->time = (/t&time/)
ncdf->lat  = (/q&lat/)
ncdf->lon  = (/q&lon/)
ncdf->q    = (/q_int/)
end

