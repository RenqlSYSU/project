;***************************************************************************
;filename:      calc_day_variance.ncl
;output_file:   UVT850500200_day_variance.nc 
;read the 27 years daily U,V,T at East Asia and then remove the annual cycle 
;get an array of (ncase,nday,nlevel,nlat,nlon),then calculate the variance in
;the dim of nday (the sample size is 27*365=9855)
;                                            by Ql_Ren
;                                           2017/12/21
;******************************************************************************
begin
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/Model_Team/F/nudgtest/mdata/UVT850500200_day_variance.nc"
var_name = (/"U","V","T"/)
case     = (/"AMIP_CTRL","TP_NUDG","DIFF"/)
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
nday = nyear*365
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
var       = new((/nvar,ncase-1,nyear,365,nlev,nlat,nlon/),float)
rmvannual = new((/nvar,ncase-1,nday,nlev,nlat,nlon/),float)
vari      = new((/nvar,ncase,nlev,nlat,nlon/),float)

;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nv = 0,nvar-1,1 
do nc = 0,ncase-2,1
    do nt = 0,nyear-1,1
        f  = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nv,nc,nt,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    end do  ;the loop of the 27 years
end do   ;the loop of the two cases 
end do   ;the loop of the 3 variables

;-------remove the annual cycle-------------
annual = dim_avg_n_Wrap(var,2) ;get (nvar,ncase-1.365,nlev,nlat,nlon)
do nt =0,nyear-1,1
   nt_s = nt*365
   nt_e = (nt+1)*365-1
   rmvannual(:,:,nt_s:nt_e,:,:,:) = var(:,:,nt,:,:,:,:) - annual
end do

vari(:,0:1,:,:,:) = dim_variance_n_Wrap(rmvannual,2)
vari(:,2,:,:,:)   = vari(:,0,:,:,:) - vari(:,1,:,:,:) ;get (nvar,ncase,nlev,nlat,nlon)
printVarSummary(vari)
prob = ftest(vari(:,0,:,:,:),nday,vari(:,1,:,:,:),nday,0)   ;(nvar,nlev,nlat,nlon)
printVarSummary(prob)
printVarSummary(annual)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The three cases are CTRL,NUDG and DIFF" + \
                    "and the three variables are U,V,T"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"var","case","lev","lat","lon"/)
dimSizes = (/nvar,ncase,nlev,nlat,nlon/)
dimUnlim = (/False,False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var","integer",(/"var"/))
filevardef(ncdf,"case","integer",(/"case"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vari",typeof(vari),(/"var","case","lev","lat","lon"/))
filevardef(ncdf,"prob",typeof(prob),(/"var","lev","lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
;vari@long_name = "variance of 27 years"
vari@units   = ""
filevarattdef(ncdf,"vari",vari)
filevarattdef(ncdf,"prob",prob)

ncdf->var    = (/1,2,3/)  ;the three variables are U,V,T
ncdf->case   = (/1,2,3/)  ;the two cases are (/"CTRL","NUDG"/)
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->vari = (/vari/)
ncdf->prob = (/prob/)
end

