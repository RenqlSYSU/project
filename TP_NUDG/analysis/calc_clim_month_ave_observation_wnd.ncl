;***************************************************************************
;filename:      calc_clim_month_ave_observation_wnd.ncl
;output_file:   Observ-Clim_month_ave_wnd.nc 
;use cpc global precip daily 1x1 data to calculate the observation climate rainfall
;which used to state the CTRL's result is right
;                                            by Ql_Ren
;                                           2018/04/18
;******************************************************************************
begin
;path    = "/users/yangsong3/renql/project/TP_NUDG"
path    = "/home/ys17-19/renql/project/TP_NUDG"
filein  = path + "/data/observation/NCEP2_daily/"  ;+year(i) + ".nc"
fileout = path + "/analysis/mdata/Observ-Clim_month_ave"
var_name = (/"uwnd","vwnd","hgt"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein + var_name(0) + "/" + var_name(0) + "." + years(0) + ".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
JJAw  = new((/nyear,nlev,nlat,nlon/),float)

do nv = 2 ,nvar-1,1 
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    f    = addfile(filein + var_name(nv) + "/" + var_name(nv) + "." + years(nt) + ".nc","r")
    time        := cd_calendar(f->time,0)
    select_time := ind(time(:,1).ge.6.and.time(:,1).le.8)
    var := short2flt(f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr}))
    JJAw(nt,:,:,:) = dim_avg_n_Wrap(var,0)
end do   ;the loop of 27 years
JJA := dim_avg_n_Wrap(JJAw,0) ;get a var(nlat,nlon)
printVarSummary(JJA)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout+var_name(nv)+".nc")
ncdf = addfile(fileout+var_name(nv)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "Climate JJA average" + var_name(nv)
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lev","lat","lon"/)
dimSizes = (/nlev ,nlat ,nlon/)
dimUnlim = (/False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,var_name(nv),typeof(JJA),(/"lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,var_name(nv),JJA)

ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->$var_name(nv)$   = (/JJA/)
end do
end

