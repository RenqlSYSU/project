;******************************************************************************
;calculate the precipitation index of zonal JJA and 
;calculate the correlation between the index and other field precipitation 
;******************************************************************************
begin
area     = "India"
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
;fileout  = "/users/yangsong3/Model_Team/F/clm-month/AMIP_C5PM." + area + "_Asia_preci_corr.nc"
fileout  = "/users/yangsong3/Model_Team/F/clm-month/AMIP_C5PM_TP_NUDG." + area + "_Asia_preci_corr.nc"
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of index
ilats = 20
ilatn = 35
ilonl = 60
ilonr = 80
;area of field variables,larger then the plan
flats = -20 ;-15
flatn = 60 ;55
flonl = 20 ;30
flonr = 160 ;150

;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
;define the dimension of the data
f    = addfile(filein+years(0)+".daily."+var_name(0)+".nc","r")
var  = f->$var_name(0)$(:,{ilats:ilatn},{ilonl:ilonr})
time = cd_calendar(var&time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)

nyear = dimsizes(years)
ntime = dimsizes(select_time)
nilat = dimsizes(var&lat)
nilon = dimsizes(var&lon)
vars  = new((/nyear,ntime,nilat,nilon/),float)
index = new((nyear),float)

;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
do t = 0,nyear-1,1
    vars(t,:,:,:) = 0
    do v = 0, dimsizes(var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
    f  = addfile(filein+years(t)+".daily."+var_name(v)+".nc","r")
    vars(t,:,:,:)  = vars(t,:,:,:) + f->$var_name(v)$(select_time,{ilats:ilatn},{ilonl:ilonr})
    end do
end do
printVarSummary(vars)
vars := dim_avg_n_Wrap(vars,1)  ;calculate the JJA average
index = wgt_areaave_Wrap(vars,1.0,1.0,0) ;calculate the area average and get the index

;====================================================================
;read the field variable, then calculate the correlation between index and area
;======================================================================
;define the dimension of the data
f    = addfile(filein+years(0)+".daily."+var_name(0)+".nc","r")
var := f->$var_name(0)$(:,{flats:flatn},{flonl:flonr})
nflat  = dimsizes(var&lat)
nflon  = dimsizes(var&lon)
fvars  = new((/nyear,ntime,nflat,nflon/),float)

do t = 0,nyear-1,1
    f_PRECC  = addfile(filein+years(t)+".daily."+var_name(0)+".nc","r")
    f_PRECL  = addfile(filein+years(t)+".daily."+var_name(1)+".nc","r")
    PRECC  = f_PRECC->$var_name(0)$(select_time,{flats:flatn},{flonl:flonr})
    PRECL  = f_PRECL->$var_name(1)$(select_time,{flats:flatn},{flonl:flonr})
    fvars(t,:,:,:) = PRECC + PRECL
end do
printVarSummary(fvars)
fvars := dim_avg_n_Wrap(fvars,1)  ;calculate the JJA average,the fvars is TLL

corr = escorc_n(index,fvars,0,0)
printVarSummary(corr)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileattdef(ncdf,fileAtt)

;Define all Coordinate
dimNames = (/"lat","lon"/)
dimSizes = (/nflat,nflon/)
dimUnlim = (/False, False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"corr",typeof(corr),getvardims(var(0,:,:)))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
copy_VarMeta(var(0,:,:),corr)
corr@long_name = "the correlation between Indian preci and " + \
                    "Asia preci (PRECC+PRECL) in JJA"
corr@units     = ""
filevarattdef(ncdf,"corr",corr)

ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr = (/corr/)
end

