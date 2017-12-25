;******************************************************************************
;calculate the four areas JJA precipitation index and 
;calculate the correlation between the index and other field precipitation 
;                                            by Ql_Ren
;                                           2017/12/02
;******************************************************************************
begin
;the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nc"
fileout  = "/users/yangsong3/Model_Team/F/clm-month/areas_telec_preci_JJA_corr.nc"
var_name = (/"PRECC","PRECL"/)
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

;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
nyear = dimsizes(years)
narea = dimsizes(area)
index = new((/2,narea,nyear/),float)

;find the time index of the JJA
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
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
        do v = 0, dimsizes(var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
            f  = addfile(filein(nf)+years(t)+".daily."+var_name(v)+".nc","r")
            vars := f->$var_name(v)$(select_time,{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
            vars := dim_avg_n_Wrap(vars,0)  ;calculate the JJA average
            index(nf,na,t) = index(nf,na,t) + wgt_areaave_Wrap(vars,1.0,1.0,0) ;calculate the area average and get the index
        end do ;the loop of precc and precl
    end do  ;the loop of the years
end do   ;the loop of the areas
end do   ;the loop of the case
printVarSummary(index)

;====================================================================
;read the field variable, then calculate the correlation between index and area
;======================================================================
;define the dimension of the data
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
var := f->$var_name(0)$(:,{flats:flatn},{flonl:flonr})
nflat  = dimsizes(var&lat)
nflon  = dimsizes(var&lon)
fvars  = new((/2,nyear,ntime,nflat,nflon/),float)

do nf = 0,1,1
do t = 0,nyear-1,1
    f_PRECC  = addfile(filein(nf)+years(t)+".daily."+var_name(0)+".nc","r")
    f_PRECL  = addfile(filein(nf)+years(t)+".daily."+var_name(1)+".nc","r")
    PRECC  = f_PRECC->$var_name(0)$(select_time,{flats:flatn},{flonl:flonr})
    PRECL  = f_PRECL->$var_name(1)$(select_time,{flats:flatn},{flonl:flonr})
    fvars(nf,t,:,:,:) = PRECC + PRECL
end do
end do
printVarSummary(fvars)
fvars := dim_avg_n_Wrap(fvars,2)  ;calculate the JJA average,the fvars is TLL

corr = new((/dimsizes(case),narea,nflat,nflon/),float)
do nf = 0,1,1
corr(nf,:,:,:) = escorc_n(index(nf,:,:),fvars(nf,:,:,:),1,0)  ;get an array of (/ncase,narea,nlat,nlon/) 
end do
corr(2,:,:,:) = corr(0,:,:,:)-corr(1,:,:,:)  ;calculate the diff
printVarSummary(corr)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The three case are CTRL,NUDG,DIFF(CTRL-NUDG)" + \
                    "the four areas are India,Tibet,Philippines,Maritime_continent"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"case","area","lat","lon"/)
dimSizes = (/dimsizes(case),dimsizes(area),nflat,nflon/)
dimUnlim = (/False, False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"case","float",(/"case"/))
filevardef(ncdf,"area","float",(/"area"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"corr",typeof(corr),(/"case","area","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
copy_VarMeta(var(0,:,:),corr(0,0,:,:))
corr@long_name = "the teleconnection of Aisa with preci (PRECC+PRECL) index in JJA"
corr@units     = ""
filevarattdef(ncdf,"corr",corr)

ncdf->case= (/1,2,3/)  ;the three case are (/"CTRL","NUDG","DIFF"/)
ncdf->area= (/1,2,3,4/) ;the four areas are (/"India","Tibet","Philippines","Maritime_continent"/)
ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr= (/corr/)
end

