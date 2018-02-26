;******************************************************************************
;filename:      calc_corr_rain_and_divg.ncl 
;output_file:   Corr_winter_rain_and_divg_vort.nc
;read ruyuan 1970-2016 rainfall and East Asia divg in winter then 
;calculation the correlation between them.
;                                            by Ql_Ren
;                                           2018/02/22
;******************************************************************************
begin
path     = "/home/ys17-19/renql/project/2017downscaling"
var_name = (/"rain","uwnd","vwnd"/)
filein1  = path + "/data/rain7016.nc"
filein2  = path + "/data/NECP_NCAR-daily/" + var_name(1:2) + "/" + var_name(1:2) + "*"
fileout  = path + "/correlation/mdata/" + "corr_winter_rain_and_divg_vort-cfd.nc" 
years    = ispan(1970,2016,1)   ;the calculate time

;area of field 
lats = 5 ;-15
latn = 60 ;55
lonl = 70 ;30
lonr = 150 ;150
lev = (/850,700,500,300/)

;==========================================================================
;read the winter(9,10,11,12,1,2) omega data,2016-1948=68,1970-1948=22 and ruyuan rainfall
;===========================================================================
files       = systemfunc("ls " + filein2(0)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
time        = cd_calendar(f2[:]->time,0)
select_time = ind(time(:,1).ge.9.or.time(:,1).le.2)
ntime       = dimsizes(select_time)
var_u      := f2[:]->$var_name(1)$(select_time,{lev},:,:)

files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v       := f2[:]->$var_name(2)$(select_time,{lev},:,:)

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

var2 = uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 3) ;calculate the divergence
var3 = uv2vr_cfd(var_u, var_v, var_v&lat, var_v&lon, 3) ;calculate the relative vorticity 
copy_VarMeta(var_u,var2)
copy_VarMeta(var_u,var3)
;var3 = var_u
;var2 = var_v
;uv2vrdvf(var_u, var_v, var3, var2) ;Computes vorticity and divergence via spherical harmonics

nlat = dimsizes(var2&lat)
nlon = dimsizes(var2&lon)
nlev = dimsizes(var2&level)

f1      = addfile(filein1,"r")
var1    = f1->$var_name(0)$(select_time)  ;get var1(time)

;=================================================================
;calculation the correlation and do the significance testing
;================================================================
corr_divg      = escorc_n(var1,var2,0,0) ;get (lev,lat,lon)
prob_divg      = rtest(corr_divg,ntime,0)
corr_vort      = escorc_n(var1,var3,0,0) ;get (lev,lat,lon)
prob_vort      = rtest(corr_vort,ntime,0)

corr_divg1day      = escorc_n(var1(1:ntime-1),var2(0:ntime-2,:,:,:),0,0) ;get (lev,lat,lon)
prob_divg1day      = rtest(corr_divg1day,ntime-1,0)
corr_vort1day      = escorc_n(var1(1:ntime-1),var3(0:ntime-2,:,:,:),0,0) ;get (lev,lat,lon)
prob_vort1day      = rtest(corr_vort1day,ntime-1,0)
printVarSummary(var2)
printVarSummary(var3)
printVarSummary(corr_divg)
printVarSummary(prob_divg)
print(ntime)
print(nlev)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the correlation between ruyuan daily rainfall(mm) and divergence as well as vorticity" + \ 
                        " in winter(9,10,11,12,1,2) from 1970 to 2016"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lev","lat" ,"lon"/)
dimSizes = (/nlev ,nlat  ,nlon /)
dimUnlim = (/False ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev",typeof(var2&level),(/"lev"/))
filevardef(ncdf,"lat",typeof(var2&lat),getvardims(var2&lat))
filevardef(ncdf,"lon",typeof(var2&lon),getvardims(var2&lon))
filevardef(ncdf,"corr_divg",typeof(corr_divg),(/"lev","lat","lon"/))
filevardef(ncdf,"corr_vort",typeof(corr_vort),(/"lev","lat","lon"/))
filevardef(ncdf,"prob_divg",typeof(prob_divg),(/"lev","lat","lon"/))
filevardef(ncdf,"prob_vort",typeof(prob_vort),(/"lev","lat","lon"/))

filevardef(ncdf,"corr_divg1day",typeof(corr_divg1day),(/"lev","lat","lon"/))
filevardef(ncdf,"corr_vort1day",typeof(corr_vort1day),(/"lev","lat","lon"/))
filevardef(ncdf,"prob_divg1day",typeof(prob_divg1day),(/"lev","lat","lon"/))
filevardef(ncdf,"prob_vort1day",typeof(prob_vort1day),(/"lev","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lev",var2&level)
filevarattdef(ncdf,"lat",var2&lat)
filevarattdef(ncdf,"lon",var2&lon)
corr_divg@long_name = "the simultaneous correlation with divergence" 
corr_divg@units     = ""
filevarattdef(ncdf,"corr_divg",corr_divg)
corr_vort@long_name = "the simultaneous correlation with vorticity" 
corr_vort@units     = ""
filevarattdef(ncdf,"corr_vort",corr_vort)
filevarattdef(ncdf,"prob_divg",prob_divg)
filevarattdef(ncdf,"prob_vort",prob_vort)

filevarattdef(ncdf,"corr_divg1day",corr_divg1day)
filevarattdef(ncdf,"corr_vort1day",corr_vort1day)
filevarattdef(ncdf,"prob_divg1day",prob_divg1day)
filevarattdef(ncdf,"prob_vort1day",prob_vort1day)

ncdf->lev = (/var2&level/) 
ncdf->lat = (/var2&lat/)
ncdf->lon = (/var2&lon/)
ncdf->corr_divg1day      = (/corr_divg1day/)
ncdf->corr_vort1day      = (/corr_vort1day/)
ncdf->prob_divg1day      = (/prob_divg1day/)
ncdf->prob_vort1day      = (/prob_vort1day/)

ncdf->corr_divg      = (/corr_divg/)
ncdf->corr_vort      = (/corr_vort/)
ncdf->prob_divg      = (/prob_divg/)
ncdf->prob_vort      = (/prob_vort/)
end
