;***************************************************************************
;filename:      calc_clim_month_ave_observation_preci.ncl
;output_file:   Observ-Clim_month_ave_preci.nc 
;use cpc global precip daily 1x1 data to calculate the observation climate rainfall
;which used to state the CTRL's result is right
;                                            by Ql_Ren
;                                           2018/04/18
;******************************************************************************
begin
;path    = "/users/yangsong3/renql/project/TP_NUDG"
path    = "/home/ys17-19/renql/project/TP_NUDG"
filein  = path + "/data/observation/gpcp_monthly_preci/precip.mon.mean.nc" ;time from 197901 to 201711, monthly data 
fileout = path + "/analysis/mdata/GPCP-Clim_month_ave_preci.nc"
var_name = (/"precip"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein,"r")
time = cd_calendar( f->time,0)
select_time = ind(time(:,0).ge.1979.and.time(:,0).le.2005.and.time(:,1).ge.6.and.time(:,1).le.8)
var  = f->$var_name$(select_time,{lats:latn},{lonl:lonr})
var := dim_avg_n_Wrap(var,0) ;get (lat,lon)
nlon = dimsizes(var&lon)
nlat = dimsizes(var&lat)
printVarSummary(var)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "Climate JJA average precipation from 1979 to 2005"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lat","lon"/)
dimSizes = (/nlat ,nlon/)
dimUnlim = (/False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"preci",typeof(var),(/"lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"preci", var)

ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->preci   = (/var/)
end

