;***************************************************************************
;filename:      calc_day_variance.ncl
;output_file:   UVT850500200_day_variance.nc 
;read the 27 years daily U,V,T at East Asia and then remove the annual cycle 
;get an array of (ncase,nday,nlevel,nlat,nlon),then calculate the variance in
;the dim of nday (the sample size is 27*365=9855)
;                                            by Ql_Ren
;                                           2017/12/21
;modified
;2018/01/08 for calculate the new case TP-NUDG-6h
;2018/01/17 for calculate one case
;******************************************************************************
begin
;the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/nudgtest/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-day_variance_preci.nc"
var_name = (/"PRECC","PRECL"/)   ;(/"U","V","T"/)
case     = (/"AMIP_CTRL","TP_NUDG24h","TP_NUDG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time
ispreci = 1

;area: east Asia
lats = -20
latn = 60
lonl = 20
lonr = 220
lev  = (/850,500,200/)

if(ispreci.eq.0) then
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nyear= dimsizes(years)
nday = nyear*365
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
var       = new((/nvar,nyear,365,nlev,nlat,nlon/),float)
rmvannual = new((/nvar,nday,nlev,nlat,nlon/),float)
vari      = new((/nvar,nlev,nlat,nlon/),float)

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
;-------remove the annual cycle-------------
    annual = dim_avg_n_Wrap(var,1) ;get (nvar,365,nlev,nlat,nlon)
    do nt =0,nyear-1,1
        nt_s = nt*365
        nt_e = (nt+1)*365-1
        rmvannual(:,nt_s:nt_e,:,:,:) = var(:,nt,:,:,:,:) - annual
    end do
    vari(:,:,:,:) = dim_variance_n_Wrap(rmvannual,1)  ;get (nvar,nlev,nlat,nlon)
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
dimNames = (/"var","lev","lat","lon"/)
dimSizes = (/nvar,nlev,nlat,nlon/)
dimUnlim = (/False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var","integer",(/"var"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vari",typeof(vari),(/"var","lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
;vari@long_name = "variance of 27 years"
vari@units   = ""
filevarattdef(ncdf,"vari",vari)

ncdf->var    = (/1,2,3/)  ;the three variables are U,V,T
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->vari = (/vari/)
end do

;-------add var_name(0) and var_name(1)-------------
else
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nyear= dimsizes(years)
nday = nyear*365
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
var       = new((/nyear,365,nlat,nlon/),float)
rmvannual = new((/nday,nlat,nlon/),float)
vari      = new((/nlat,nlon/),float)

do nc = 0,ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
    var = 0
    do nt = 0,nyear-1,1
    do nv = 0,nvar-1,1 
        f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do  ;the loop of the PRECC and PRECL
    end do   ;the loop of 27 year
;-------remove the annual cycle-------------
    var = var*1000*3600*24
    annual = dim_avg_n_Wrap(var,0) ;get (365,nlat,nlon)
    do nt =0,nyear-1,1
        nt_s = nt*365
        nt_e = (nt+1)*365-1
        rmvannual(nt_s:nt_e,:,:) = var(nt,:,:,:) - annual
    end do
    vari(:,:) = dim_variance_n_Wrap(rmvannual,0)  ;get (nlev,nlat,nlon)
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
fileAtt@discription   = "preci"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lat","lon"/)
dimSizes = (/nlat,nlon/)
dimUnlim = (/False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vari",typeof(vari),(/"lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
vari@long_name = "daily variance of preci"
vari@units   = ""
filevarattdef(ncdf,"vari",vari)

ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->vari = (/vari/)
end do

end if
end

