;***************************************************************************
;read the 27 years U at East Asia and then calculate the month average to 
;get an array of (ncase,nyear,nmonth,nlat,nlon),then calculate the variance in
;the dim of nyear (the sample size is 27)
;                                            by Ql_Ren
;                                           2017/12/07
;******************************************************************************
begin
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/Model_Team/F/nudgtest/mdata/Month_ave_year_variance.nc"
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
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
nlev  = dimsizes(lev)
var       = new((/nvar,ncase,nyear,365,nlev,nlat,nlon/),float)
var_month = new((/nvar,ncase,nyear,12,nlev,nlat,nlon/),float)
vari      = new((/nvar,ncase,12,nlev,nlat,nlon/),float)

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

time = cd_calendar(f->time,0)
do nm = 0,11,1
    select_time := ind(time(:,1).eq.nm+1)
    var_month(:,0:1,:,nm,:,:,:) = dim_avg_n_Wrap(var(:,0:1,:,select_time,:,:,:),3)
end do

vari(:,0:1,:,:,:,:) = dim_variance_n_Wrap(var_month(:,0:1,:,:,:,:,:),2)
vari(:,2,:,:,:,:)   = vari(:,0,:,:,:,:) - vari(:,1,:,:,:,:) ;(nvar,ncase,12,nlev,nlat,nlon)
printVarSummary(vari)
prob = ftest(vari(:,0,:,:,:,:),nyear,vari(:,1,:,:,:,:),nyear,0)   ;(nvar,12,nlev,nlat,nlon)
printVarSummary(prob)

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
dimNames = (/"var","case","month","lev","lat","lon"/)
dimSizes = (/nvar,ncase,12,nlev,nlat,nlon/)
dimUnlim = (/False, False,False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var","integer",(/"var"/))
filevardef(ncdf,"case","integer",(/"case"/))
filevardef(ncdf,"month","integer",(/"month"/))
filevardef(ncdf,"lev","integer",(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vari",typeof(vari),(/"var","case","month","lev","lat","lon"/))
filevardef(ncdf,"prob",typeof(prob),(/"var","month","lev","lat","lon"/))

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
ncdf->month  = ispan(1,12,1)
ncdf->lev = (/lev/)
ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->vari = (/vari/)
ncdf->prob = (/prob/)
end

