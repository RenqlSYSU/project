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
filein  = path + "/data/observation/cpc_daily_preci/precip."  ;+year(i) + ".nc"
fileout = path + "/analysis/mdata/Observ-Clim_month_ave_preci.nc"
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
f    = addfile(filein+years(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
JJApreci  = new((/nyear,nlat,nlon/),float)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f    = addfile(filein+years(nt)+".nc","r")
    time        := cd_calendar(f->time,0)
    select_time := ind(time(:,1).ge.6.and.time(:,1).le.8)
    var := f->$var_name$(select_time,{lats:latn},{lonl:lonr})
    JJApreci(nt,:,:) = dim_avg_n_Wrap(var,0)
end do   ;the loop of 27 years
JJApreci := dim_avg_n_Wrap(JJApreci,0) ;get a var(nlat,nlon)
printVarSummary(JJApreci)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "Climate JJA average precipation"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lat","lon"/)
dimSizes = (/nlat ,nlon/)
dimUnlim = (/False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"preci",typeof(JJApreci),(/"lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"preci",JJApreci)

ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->preci   = (/JJApreci/)
end

