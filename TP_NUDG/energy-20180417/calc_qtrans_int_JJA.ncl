;***************************************************************************
;filename:      calc_qtran_int_JJA.ncl
;use the initial data to calculate the integrated qu and qv then write data file (a year is a file)
;The ncl functions have been used are: dpres_hybrid_ccm,dim_sum_n_Wrap,dim_avg_n_Wrap
;                                            by Ql_Ren
;                                           2018/04/17
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = (/"/users/yangsong3/L_Zealot/F/AMIP_C5PM/exe/AMIP_C5PM.cam.h1.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/pre_data_24h/AMIP_C5PM_TP_NUDG.cam.h1.",\
            "/users/yangsong3/L_Zealot/F/AMIP_C5PM_TP_NUDG/exe/AMIP_C5PM_TP_NUDG.cam.h1."/)
fileout = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
case    = (/"CTRL","NUDG24h","NUDG6h"/)
years   = ispan(1979,2005,1)   ;the calculate time
ncase = dimsizes(case)
nyear = dimsizes(year)

files1 = systemfunc("ls " + filein(0) + years(0) + "*" + ".nc")
f1     = addfile(files1(0),"r")
hyai   = f1->hyai  ;hyai(ilev),ilev = lev+1,hybrid A coefficient at layer interfaces
hybi   = f1->hybi
hyam   = f1->hyam  ;hybm(lev), hybrid A coefficient at layer midpoints
hybm   = f1->hybm
p0     = 100000   ;a scalar ,double,reference pressure,unit is Pa
printVarSummary(hyai)
printVarSummary(hyam)

do nf = 0,ncase-1,1 
do nt = 0,nyear-1,1
;==========================================================
;read the data
;============================================================
files  = systemfunc("ls " + filein(nf) + years(nt) + "*" + ".nc")
f      = addfiles(files,"r")
ListSetType(f,"cat")
u   = f[:]->U    ;(:,:,{lats:latn},{lonl:lonr})
v   = f[:]->V    ;(:,:,{lats:latn},{lonl:lonr})
q   = f[:]->Q    ;(:,:,{lats:latn},{lonl:lonr})
ps  = f[:]->PS   ;(:,{lats:latn},{lonl:lonr})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
printVarSummary(ps)
printVarSummary(u)
printVarSummary(q)
q = q*1000 ;convert unit from kg/kg to g/kg

;calculates the pressure differences of a hybrid coordinate system.
dpi  = dpres_hybrid_ccm(ps,p0,hyai,hybi)  ;get dpi (time,lev,lat,lon)
printVarSummary(dpi)

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

