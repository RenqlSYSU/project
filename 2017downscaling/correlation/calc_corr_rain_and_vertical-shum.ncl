;******************************************************************************
;filename:      calc_corr_rain_and_vertical-shum.ncl 
;output_file:   corr_winter_rain_and_vertical-shum.nc
;read ruyuan 1970-2016 rainfall and global shum and shum flux divergence of 
;vertical sum  in winter then calculation the correlation between them.
;                                            by Ql_Ren
;                                           2018/02/24
;******************************************************************************
begin
path     = "/home/ys17-19/renql/project/2017downscaling"
var_name = (/"rain","uwnd","vwnd","shum","pres.sfc"/)
filein1  = path + "/data/rain7016.nc"
filein2  = path + "/data/NECP_NCAR-daily/" + var_name(1:4) + "/" + var_name(1:4) + "*"
fileout  = path + "/correlation/mdata/" + "corr_winter_rain_and_vertical-shum-cfd.nc" 
years    = ispan(1970,2016,1)   ;the calculate time
g        = 9.8 ;m/(s*s)

;area of field 
lats = 5 ;-15
latn = 60 ;55
lonl = 70 ;30
lonr = 150 ;150
;lev = (/850,700,500,300/)

;==========================================================================
;read the winter(9,10,11,12,1,2) omega data,2016-1948=68,1970-1948=22 and ruyuan rainfall
;===========================================================================
files       = systemfunc("ls " + filein2(0)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
time        = cd_calendar(f2[:]->time,0)
select_time = ind(time(:,1).ge.9.or.time(:,1).le.2)
ntime       = dimsizes(select_time)
var_u      := f2[:]->$var_name(1)$(select_time,0:7,:,:)

files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v       := f2[:]->$var_name(2)$(select_time,0:7,:,:)

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var        := f2[:]->$var_name(3)$(select_time,:,:,:)
var_shum    = var*1000      ;convert unit from kg/kg to g/kg
copy_VarMeta(var,var_shum)

files       = systemfunc("ls " + filein2(3)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var        := f2[:]->pres(select_time,:,:)
pres_sfc    = var/100       ;convert unit from Pa to hPa
copy_VarMeta(var,pres_sfc)

var_flux    = wgt_vert_avg_beta ( var_shum&level, var_shum, pres_sfc, 0, 0 ) / g 
;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum

if(any(ismissing(var_u))) then
    print("uwnd contains some missing values. Beware.")
else
    print("uwnd don't have missing values, so can use spherical harmonics")
end if
if(any(ismissing(var_v))) then
    print("vwnd contains some missing values. Beware.")
else
    print("vwnd don't have missing values, so can use spherical harmonics")
end if
if(any(ismissing(var_shum))) then
    print("shum contains some missing values. Beware.")
else
    print("shum don't have missing values, so can use spherical harmonics")
end if

var_u    = var_u*var_shum/g
var_v    = var_v*var_shum/g
var_divg = uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 3) ;calculate the divergence
;var_divg = uv2dvF_Wrap( var_u, var_v )

var_divg_flux    = wgt_vert_avg_beta ( var_shum&level, var_divg, pres_sfc, 0, 0 )  
;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum

printVarSummary(var_flux)
printVarSummary(var_divg_flux)

nlat = dimsizes(var_u&lat)
nlon = dimsizes(var_u&lon)

f1      = addfile(filein1,"r")
var1    = f1->$var_name(0)$(select_time)  ;get var1(time)

;=================================================================
;calculation the correlation and do the significance testing
;================================================================
corr_flux      = escorc_n(var1,var_flux,0,0) ;get (lev,lat,lon)
prob_flux      = rtest(corr_flux,ntime,0)
corr_divg_flux      = escorc_n(var1,var_divg_flux,0,0) ;get (lev,lat,lon)
prob_divg_flux      = rtest(corr_divg_flux,ntime,0)


corr_flux1day      = escorc_n(var1(1:ntime-1),var_flux(0:ntime-2,:,:),0,0) ;get (lev,lat,lon)
prob_flux1day      = rtest(corr_flux1day,ntime-1,0)
corr_divg_flux1day      = escorc_n(var1(1:ntime-1),var_divg_flux(0:ntime-2,:,:),0,0) ;get (lev,lat,lon)
prob_divg_flux1day      = rtest(corr_divg_flux1day,ntime-1,0)

corr_flux@units      = "g/(m*m)"
corr_flux@long_name  = "the correlation between ruyuan rain and The layer of water vapor from 1000mb to 300mb"
corr_divg_flux@units      = "g/(s*m*m)"
corr_divg_flux@long_name  = "the correlation between ruyuan rain and The layer of water vapor divergence from 1000mb to 300mb"

printVarSummary(corr_flux)
printVarSummary(prob_flux)
print(ntime)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the correlation between ruyuan daily rainfall(mm) and the layer of water vapor as well as water vapor divergence" + \ 
                        " in winter(9,10,11,12,1,2) from 1970 to 2016"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lat" ,"lon"/)
dimSizes = (/nlat  ,nlon /)
dimUnlim = (/False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))

filevardef(ncdf,"corr_flux",typeof(corr_flux),(/"lat","lon"/))
filevardef(ncdf,"corr_divg_flux",typeof(corr_divg_flux),(/"lat","lon"/))
filevardef(ncdf,"prob_flux",typeof(prob_flux),(/"lat","lon"/))
filevardef(ncdf,"prob_divg_flux",typeof(prob_divg_flux),(/"lat","lon"/))

filevardef(ncdf,"corr_flux1day",typeof(corr_flux1day),(/"lat","lon"/))
filevardef(ncdf,"corr_divg_flux1day",typeof(corr_divg_flux1day),(/"lat","lon"/))
filevardef(ncdf,"prob_flux1day",typeof(prob_flux1day),(/"lat","lon"/))
filevardef(ncdf,"prob_divg_flux1day",typeof(prob_divg_flux1day),(/"lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)

filevarattdef(ncdf,"corr_flux",corr_flux)
filevarattdef(ncdf,"corr_divg_flux",corr_divg_flux)
filevarattdef(ncdf,"prob_flux",prob_flux)
filevarattdef(ncdf,"prob_divg_flux",prob_divg_flux)

filevarattdef(ncdf,"corr_flux1day",corr_flux1day)
filevarattdef(ncdf,"corr_divg_flux1day",corr_divg_flux1day)
filevarattdef(ncdf,"prob_flux1day",prob_flux1day)
filevarattdef(ncdf,"prob_divg_flux1day",prob_divg_flux1day)

ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr_flux      = (/corr_flux/)
ncdf->corr_divg_flux      = (/corr_divg_flux/)
ncdf->prob_flux      = (/prob_flux/)
ncdf->prob_divg_flux      = (/prob_divg_flux/)

ncdf->corr_flux1day      = (/corr_flux1day/)
ncdf->corr_divg_flux1day      = (/corr_divg_flux1day/)
ncdf->prob_flux1day      = (/prob_flux1day/)
ncdf->prob_divg_flux1day      = (/prob_divg_flux1day/)
end
