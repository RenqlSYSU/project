;******************************************************************************
;filename:      calc_Asia_JJAu_EOF.ncl
;output_file:   Asia_JJAu_covEOFn3.nc
;calculate the Asia JJA average U and 
;use the function of eofunc_n_Wrap, eofunc_ts_n_Wrap to 
;calculate the EOF pattern and time serises.
;                                            by Ql_Ren
;                                           2017/12/19
;******************************************************************************
begin
;------the path of the file
filein   = (/"/users/yangsong3/Model_Team/F/AMIP-CTRL/pro/AMIP_C5PM.cam.h1.YEAR.",\ ;+year+".daily.nc"
             "/users/yangsong3/Model_Team/F/TP-NUDG/pro/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/) ;+year+".daily.nfileout  = "/users/yangsong3/Model_Team/F/clm-month/mdata/Areas_JJApreci_EOFn3.nc"
fileout  = "/users/yangsong3/Model_Team/F/clm-month/mdata/Asia_JJAv_covEOFn3.nc"
years    = ispan(1979,2005,1)   ;the calculate time
case     = (/"CTRL","NUDG"/)
var_name = (/"V"/)
;-----area of east asia 
lats = 30
latn = 50
lonl = 30
lonr = 130
lev  = 200
neval = 4
nyear = dimsizes(years)
ncase = dimsizes(case)

;------find the time index of the JJA and define the dimension of the varreas_telec_preci_JJA_corr.nc
f           = addfile(filein(0)+years(0)+".daily."+var_name+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
var  := f->$var_name$(0,0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(var&lat)
nlon  = dimsizes(var&lon)
vars  = new((/nyear,ntime,nlat,nlon/),float)
eof   = new((/ncase,neval,nlat,nlon/),float)
eof_ts= new((/ncase,neval,nyear/),float)
pcvar = new((/ncase,neval/),float)
eval  = new((/ncase,neval/),float)

;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
opt = True
opt@jopt = 0
do nf = 0,ncase-1,1 
    do t = 0,nyear-1,1
        f              = addfile(filein(nf)+years(t)+".daily."+var_name+".nc","r")
        vars(t,:,:,:) = f->$var_name$(select_time,{lev},{lats:latn},{lonl:lonr})
    end do  ;the loop of the years
    avg_var     = dim_avg_n_Wrap(vars,1)                     ;get the (nyear,nlat,nlon)
    eof_temp    = eofunc_n_Wrap(avg_var,neval,opt,0)         ;get the (neval,nlat,nlon)
    eof_ts_temp = eofunc_ts_n_Wrap(avg_var,eof_temp,False,0) ;get the (neval,nyear)
    sig         = eofunc_north(eof_temp@pcvar,nyear,True)
 eof(nf,:,:,:)  = eof_temp                                 ;get the (ncase,neval,nlat,nlon) 
 eof_ts(nf,:,:) = dim_standardize_n_Wrap(eof_ts_temp,1,1)  ;get the (ncase,neval,nyear)
 pcvar(nf,:)    = eof_temp@pcvar
 eval(nf,:)     = eof_temp@eval
end do   ;the loop of the case
printVarSummary(eof)
printVarSummary(eof_ts)
printVarSummary(pcvar)
printVarSummary(eval)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "The two cases are CTRL,NUDG" + \
                    " the stored data is EOFs and eof_ts of JJA U."
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"neval","case","year","lat","lon"/)
dimSizes = (/neval,ncase,nyear,nlat,nlon/)
dimUnlim = (/False, False,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"case","integer",(/"case"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"eof",typeof(eof),(/"case","neval","lat","lon"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"case","neval","year"/))
filevardef(ncdf,"pcvar",typeof(pcvar),(/"case","neval"/))
filevardef(ncdf,"eval",typeof(eval),(/"case","neval"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
;corr@long_name = "the regression patterns for the 850hPa wind calculated by JJA preci index"
;corr@units     = ""
filevarattdef(ncdf,"eof",eof)
filevarattdef(ncdf,"eof_ts",eof_ts)

ncdf->neval  = (/1,2,3,4/)  
ncdf->case   = (/1,2/)  ;the three case are (/"CTRL","NUDG"/)
ncdf->year   = (/years/) 
ncdf->lat    = (/var&lat/)
ncdf->lon    = (/var&lon/)
ncdf->eof    = (/eof/)
ncdf->eof_ts = (/eof_ts/)
ncdf->pcvar  = (/pcvar/)
ncdf->eval   = (/eval/)
end

