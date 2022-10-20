;******************************************************************************
;filename:      calc_corr_rain_and_omega.ncl 
;output_file:   Corr_rain_and_omega.nc
;read ruyuan 1970-2016 rainfall and East Asia omega in winter then 
;calculation the correlation between them.
;                                            by Ql_Ren
;                                           2018/02/18
;******************************************************************************
begin
path     = "/home/ys17-19/renql/project/2017downscaling"
var_name = (/"rain","omega"/)
filein1  = path + "/data/rain7016.nc"
filein2  = path + "/data/NECP_NCAR-daily/" + var_name(1) + "/" + var_name(1) + "*"
;fileout  = path + "/correlation/mdata/" + "corr_winter_rain_and_omega.nc" 
fileout  = path + "/correlation/mdata/" + "corr_rain_and_omega.nc" 
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
files       = systemfunc("ls " + filein2) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
time        = cd_calendar(f2[:]->time,0)
;select_time = ind(time(:,1).ge.9.or.time(:,1).le.2)    ;winter
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)    ;summer
ntime       = dimsizes(select_time)
var2       := f2[:]->$var_name(1)$(select_time,{lev},:,:)
nlat = dimsizes(var2&lat)
nlon = dimsizes(var2&lon)
nlev = dimsizes(var2&level)

f1      = addfile(filein1,"r")
var1    = f1->$var_name(0)$(select_time)  ;get var1(time)

;=================================================================
;calculation the correlation and do the significance testing
;================================================================
corr      = escorc_n(var1,var2,0,0) ;get (lev,lat,lon)
corr1day  = escorc_n(var1(1:ntime-1),var2(0:ntime-2,:,:,:),0,0) ;get (lev,lat,lon)
corr2day  = escorc_n(var1(2:ntime-1),var2(0:ntime-3,:,:,:),0,0) ;get (lev,lat,lon)
prob      = rtest(corr,ntime,0)
prob1day  = rtest(corr1day,ntime-1,0)
prob2day  = rtest(corr2day,ntime-2,0)
printVarSummary(corr)
printVarSummary(prob)
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
fileAtt@discription   = "the correlation between ruyuan daily rainfall(mm) and east asia omega(pa/s)" + \ 
                        "from 1970 to 2016"
                        ;" in summer (3,4,5,6,7,8) from 1970 to 2016"
                        ;" in winter(9,10,11,12,1,2) from 1970 to 2016"
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
filevardef(ncdf,"corr",typeof(corr),(/"lev","lat","lon"/))
filevardef(ncdf,"corr1day",typeof(corr1day),(/"lev","lat","lon"/))
filevardef(ncdf,"corr2day",typeof(corr2day),(/"lev","lat","lon"/))
filevardef(ncdf,"prob",typeof(prob),(/"lev","lat","lon"/))
filevardef(ncdf,"prob1day",typeof(prob1day),(/"lev","lat","lon"/))
filevardef(ncdf,"prob2day",typeof(prob2day),(/"lev","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lev",var2&level)
filevarattdef(ncdf,"lat",var2&lat)
filevarattdef(ncdf,"lon",var2&lon)
corr@long_name = "the simultaneous correlation" 
corr@units     = ""
filevarattdef(ncdf,"corr",corr)
corr1day@long_name = "the one day delay correlation" 
corr1day@units     = ""
filevarattdef(ncdf,"corr1day",corr1day)
corr2day@long_name = "the two day delay correlation" 
corr2day@units     = ""
filevarattdef(ncdf,"corr2day",corr2day)
filevarattdef(ncdf,"prob",prob)
filevarattdef(ncdf,"prob1day",prob1day)
filevarattdef(ncdf,"prob2day",prob2day)

ncdf->lev = (/var2&level/) 
ncdf->lat = (/var2&lat/)
ncdf->lon = (/var2&lon/)
ncdf->corr      = (/corr/)
ncdf->corr1day  = (/corr1day/)
ncdf->corr2day  = (/corr2day/)
ncdf->prob      = (/prob/)
ncdf->prob1day  = (/prob1day/)
ncdf->prob2day  = (/prob2day/)
end
