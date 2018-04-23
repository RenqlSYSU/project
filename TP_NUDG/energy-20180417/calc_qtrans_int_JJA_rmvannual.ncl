;***************************************************************************
;filename:      calc_qtran_int_JJA.ncl
;use the initial data to calculate the integrated qu and qv then write data file (a year is a file)
;The ncl functions have been used are: dpres_hybrid_ccm,dim_sum_n_Wrap,dim_avg_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_month_ave_qtran"
var_name= (/"U","V","Q","PS"/)
case    = (/"CTRL","NUDG24h","NUDG6h"/)
years   = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/1000,925,850,700,500/)

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

print("convert unit from kg/kg to g/kg")
q = q*1000 

;calculates the pressure differences of a hybrid coordinate system.
dp  = dpres_plevel_Wrap(lev*100,ps,50000,0)  ;unit is Pa, get dp (time,lev,lat,lon)
printVarSummary(dp)



qu = q*u
qv = q*v
qdivg = uv2dv_cfd ( qv, qu, u&lat, u&lon, 3)

qu := dim_avg_n_Wrap(qu*dpi,1)/9.8 ;get qu(time,lat,lon)
qv := dim_avg_n_Wrap(qv*dpi,1)/9.8 ;get qv(time,lat,lon)
copy_VarMeta(q(0,0,:,:),qu(0,:,:))
copy_VarMeta(q(0,0,:,:),qv(0,:,:))
qu@long_name = "Integrated zonal moisture transport"
qv@long_name = "Integrated meridional moisture transport"
qdivg@long_name = "Integrated divergence of moisture flux" 
qu@units     = "g/(s*m)" 
qv@units     = "g/(s*m)" 
qdivg@units     = "g/(s*m*m)" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+years(nt) + ".daily." + "qtran.nc")
ncdf = addfile(fileout(nf)+years(nt) + ".daily." + "qtran.nc","c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "the integration is based on 30 hybrid levels "
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lat" ,"lon"/)
dimSizes = (/ntime ,nlat  ,nlon /)
dimUnlim = (/True  ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(uwnd&time),(/"time"/))
filevardef(ncdf,"lat" ,typeof(uwnd&lat),getvardims(uwnd&lat))
filevardef(ncdf,"lon" ,typeof(uwnd&lon),getvardims(uwnd&lon))
filevardef(ncdf,"kpkm",typeof(kpkm),(/"time","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",uwnd&lat)
filevarattdef(ncdf,"lat" ,uwnd&lat)
filevarattdef(ncdf,"lon" ,uwnd&lon)
kpkm@long_name = "vertical integrated and ten grids average, >0 mean Kp convert to Km"
kpkm@units     = "J"
filevarattdef(ncdf,"kpkm",kpkm)

ncdf->time = (/uwnd&time/)
ncdf->lat  = (/uwnd&lat/)
ncdf->lon  = (/uwnd&lon/)
ncdf->kpkm = (/kpkm/)
end do
end

