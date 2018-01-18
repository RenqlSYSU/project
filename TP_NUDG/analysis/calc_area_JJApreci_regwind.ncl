;******************************************************************************
;calculate the four areas JJA precipitation index and 
;use they to regress the 850hPa winds (use the function "regCoef_n(x,y,dim_x,dim_y)")
;the nyears is the number of samples
;                                            by Ql_Ren
;                                           2017/12/02
;******************************************************************************

begin
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/Model_Team/F/clm-month/Areas_preci_JJA_reg850wind.nc"
index_var_name = (/"PRECC","PRECL"/)   ;the index is calculated by the precipitation
var_name =(/"U","V"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"India","Tibet","Philippines","Maritime_continent"/)
case  = (/"CTRL","NUDG","DIFF"/)
ilats = (/20,25,5,-10/)
ilatn = (/35,35,20,7.5/)
ilonl = (/60,90,120,85/)
ilonr = (/80,103,140,120/)

;area of field variables,larger then the plan
flats = -20 ;-15
flatn = 60 ;55
flonl = 20 ;30
flonr = 160 ;150
flev  = 850

;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
nyear = dimsizes(years)
narea = dimsizes(area)
index = new((/2,narea,nyear/),float)

;find the time index of the JJA
f    = addfile(filein(0)+years(0)+".daily."+index_var_name(0)+".nc","r")
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime = dimsizes(select_time)

;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
do nf = 0,1,1 
do na=0,narea-1,1
    do t = 0,nyear-1,1
        index(nf,na,t) = 0
        do v = 0, dimsizes(index_var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
            f  = addfile(filein(nf)+years(t)+".daily."+index_var_name(v)+".nc","r")
            vars := f->$index_var_name(v)$(select_time,{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
            vars := dim_avg_n_Wrap(vars,0)  ;calculate the JJA average
            index(nf,na,t) = index(nf,na,t) + wgt_areaave_Wrap(vars,1.0,1.0,0) ;calculate the area average and get the index
        end do ;the loop of precc and precl
    end do  ;the loop of the years
end do   ;the loop of the areas
end do   ;the loop of the case
printVarSummary(index)
index = dim_standardize_n_Wrap(index,1,2)

;====================================================================
;read the field variable, then calculate the correlation between index and area
;======================================================================
;define the dimension of the data
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
var := f->$var_name(0)$(0,{flev},{flats:flatn},{flonl:flonr})  ;get the var(lat,lon)
nflat  = dimsizes(var&lat)
nflon  = dimsizes(var&lon)
nvars  = dimsizes(var_name)
fvars  = new((/2,nyear,nvars,ntime,nflat,nflon/),float)  ;the 2 means AMIP-CTRL, TP-NUDG

do nf = 0,1,1
do t = 0,nyear-1,1
    do v = 0,nvars-1,1
        f    := addfile(filein(nf)+years(t)+".daily."+var_name(v)+".nc","r")
        fvars(nf,t,v,:,:,:)= f->$var_name(v)$(select_time,{flev},{flats:flatn},{flonl:flonr})
    end do
end do
end do
printVarSummary(fvars)
fvars := dim_avg_n_Wrap(fvars,3)  ;calculate the JJA average,the fvars is TLL

ncase = dimsizes(case)
corr = new((/ncase,narea,nvars,nflat,nflon/),float)
tval = new((/ncase,narea,nvars,nflat,nflon/),float)
yintercept = new((/ncase,narea,nvars,nflat,nflon/),float)
rstd = new((/ncase,narea,nvars,nflat,nflon/),float)
do nf = 0,1,1
corr(nf,:,:,:,:) = regCoef_n(index(nf,:,:),fvars(nf,:,:,:,:),1,0)  ;get an array of (/ncase,narea,nvars,nlat,nlon/) 
tval(nf,:,:,:,:) = onedtond(corr@tval,(/narea,nvars,nflat,nflon/))
yintercept(nf,:,:,:,:) = onedtond(corr@yintercept,(/narea,nvars,nflat,nflon/))
rstd(nf,:,:,:,:) = onedtond(corr@rstd,(/narea,nvars,nflat,nflon/))
end do
corr(2,:,:,:,:) = corr(0,:,:,:,:)-corr(1,:,:,:,:)  ;calculate the diff
tval(2,:,:,:,:) = tval(0,:,:,:,:)-tval(1,:,:,:,:)  ;calculate the diff
yintercept(2,:,:,:,:) = yintercept(0,:,:,:,:)-yintercept(1,:,:,:,:)  ;calculate the diff
rstd(2,:,:,:,:) = rstd(0,:,:,:,:)-rstd(1,:,:,:,:)  ;calculate the diff
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
fileAtt@discription   = "The three case are CTRL,NUDG,DIFF(CTRL-NUDG)" + \
                    "the four areas are India,Tibet,Philippines,Maritime_continent" + \
                    " the two vars is U and V"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"case","area","vars","lat","lon"/)
dimSizes = (/dimsizes(case),dimsizes(area),nvars,nflat,nflon/)
dimUnlim = (/False, False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"case","float",(/"case"/))
filevardef(ncdf,"area","float",(/"area"/))
filevardef(ncdf,"vars","float",(/"vars"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"corr",typeof(corr),(/"case","area","vars","lat","lon"/))
filevardef(ncdf,"tval",typeof(corr),(/"case","area","vars","lat","lon"/))
filevardef(ncdf,"yintercept",typeof(corr),(/"case","area","vars","lat","lon"/))
filevardef(ncdf,"rstd",typeof(corr),(/"case","area","vars","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
corr@long_name = "the regression patterns for the 850hPa wind calculated by JJA preci index"
corr@units     = ""
delete(corr@nptxy)
delete(corr@tval)
delete(corr@yintercept)
delete(corr@rstd)
filevarattdef(ncdf,"corr",corr)

ncdf->case= (/1,2,3/)  ;the three case are (/"CTRL","NUDG","DIFF"/)
ncdf->area= (/1,2,3,4/) ;the four areas are (/"India","Tibet","Philippines","Maritime_continent"/)
ncdf->vars= (/1,2/) ;the two variables is U, V
ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr= (/corr/)
ncdf->tval= (/tval/)
ncdf->yintercept= (/yintercept/)
ncdf->rstd = (/rstd/)
end

