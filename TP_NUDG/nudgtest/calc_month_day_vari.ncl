;***************************************************************************
;filename:      calc_month_day_variance.ncl
;output_file:   Month_day_variance-UVT.nc 
;read the 27 years daily U,V,T at East Asia and then
;get an array of (nvar,nyear,nday,nlevel,nlat,nlon),then calculate the variance in
;the dim of nday to get (nvar,nyear,nmonth,nlevel,nlat,nlon),then calculate
;the average in the dim of nyear.
;                                            by Ql_Ren
;                                           2018/01/17
;******************************************************************************
begin
;the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/nudgtest/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Month_day_variance_UVT.nc"
var_name = (/"U","V","T"/)
case     = (/"AMIP_CTRL","TP_NUDG24h","TP_NUDG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/850,500,200/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nyear= dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
var       = new((/nvar,nyear,365,nlev,nlat,nlon/),float)
vari_temp = new((/nvar,nyear,12,nlev,nlat,nlon/),float)
vari      = new((/nvar,12,nlev,nlat,nlon/),float)

do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nv = 0,nvar-1,1 
    do nt = 0,nyear-1,1
        f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nv,nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    end do  ;the loop of the 27 years
end do   ;the loop of the 3 variables

;-------calc the daily variance in 12 month-------------
time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    vari_temp(:,:,nm,:,:,:) = dim_variance_n_Wrap(var(:,:,select_time,:,:,:),2) ;get (nvar,nyear,nmonth,nlev,nlat,nlon)
end do

vari = dim_avg_n_Wrap(vari_temp,1)  ;get (nvar,nmonth,nlev,nlat,nlon)
printVarSummary(vari)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the three variables are U,V,T"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"var","month","lev","lat","lon"/)
dimSizes = (/nvar,12,nlev,nlat,nlon/)
dimUnlim = (/False,False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var","integer",(/"var"/))
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vari",typeof(vari),(/"var","month","lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
;vari@long_name = "variance of 27 years"
vari@units   = ""
filevarattdef(ncdf,"vari",vari)

ncdf->var    = (/1,2,3/)  ;the three variables are U,V,T
ncdf->month  = ispan(1,12,1)
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->vari = (/vari/)
end do
end

