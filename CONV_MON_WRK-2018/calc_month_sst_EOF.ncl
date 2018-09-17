;******************************************************************************
;filename:      calc_month_sst_EOF.ncl
;output_file:   var_name(nv) + area + _EOFn3.nc
;----------------------------------------------------------------------
; Concepts illustrated:
;   - Read the monthly TS(surface temperature,K) over last 30 years 
;   - use eofunc_n_Wrap,eofunc_ts_n_Wrap to calc EOF pattern and time serises
;   - store the data
;                                            by Ql_Ren
;                                           2018/09/17
;******************************************************************************
begin
unit     = (/"~S~o~N~C","hPa","W/m2","W/m2"/)    ;this can be used as right string 
var_name = (/"TS","PSL","FLNT","FSNT"/)
case     = (/"MX_AQUAP_SCYC","SF_AQUAP_SCYC","NO_TOPO_SCYC","REAL_WORLD_SCYC"/);" NG24H",THIS CAN BE USED AS LEFT STRING
pre_case = (/"(a","(b","(c","(d"/)
nv = 0

;-----EOF area-------------- 
lats = -30
latn = 30
lonl = 120
lonr = 260
;lev  = 200
neval = 4

;----path of the file and output file-----------------------------
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018" 
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fileout  = path + "/mdata/" + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E_corrEOF" + neval + ".nc"

;----Defining array dimensions-----------------
f    = addfile(filename(0) + var_name(nv) + ".1-100.nc","r")
vars = f->$var_name(nv)$(0,{lats:latn},{lonl:lonr})
nlat   = dimsizes(vars&lat) 
nlon   = dimsizes(vars&lon)
ncase  = dimsizes(case)
nyear  = 30
nmonth = 12
ntime  = nyear*nmonth
nd_s   = 1200-nyear*nmonth
nd_e   = nd_s+nyear*nmonth-1
var   = new((/ncase,ntime,nlat,nlon/),float)
eof   = new((/ncase,neval,nlat,nlon/),float)
eof_ts= new((/ncase,neval,ntime/),float)
pcvar = new((/ncase,neval/),float)
eval  = new((/ncase,neval/),float)

;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
opt = True
opt@jopt = 1 ;0 mean use the covariance matrix to compute EOFs (default), 1 mean use the correlation matrix
do nc = 0,ncase-1,1
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,:,:,:) = f->$var_name(nv)$(nd_s:nd_e,{lats:latn},{lonl:lonr})
    var(nc,:,:,:) = var(nc,:,:,:) - 273.15  ;convert unit from K to C

;remove the annual cycle
    do nm = 0, nmonth-1,1
        ave_var = dim_avg_n_Wrap(var(nc,nm:359:12,:,:),0)
        printVarSummary(ave_var)
        printVarSummary(var(nc,nm:359:12,:,:))
        var(nc,nm:359:12,:,:) = var(nc,nm:359:12,:,:) - conform(var(nc,nm:359:12,:,:),ave_var,(/1,2/))  
    end do

    eof_temp    = eofunc_n_Wrap(var(nc,:,:,:),neval,opt,0)         ;get the (neval,nlat,nlon)
    eof_ts_temp = eofunc_ts_n_Wrap(var(nc,:,:,:),eof_temp,False,0) ;get the (neval,nyear)
    sig         = eofunc_north(eof_temp@pcvar,ntime,True) ; will print the result
    
    eof(nc,:,:,:)  = eof_temp                                 ;get the (ncase,neval,nlat,nlon) 
    eof_ts(nc,:,:) = eof_ts_temp  ;get the (ncase,neval,nyear)
    pcvar(nc,:)    = eof_temp@pcvar     ;the percent variance associated with each eigenvalue.
    eval(nc,:)     = eof_temp@eval      ;eigenvalues

;test Whether the standard deviation and the eigenvalue are equal
    print(dim_stddev_n(eof_ts_temp,1)*dim_stddev_n(eof_ts_temp,1))
    print(eof_temp@eval)
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
fileAtt@discription   = "The four case, and optEOF@jopt=" + opt@jopt + \
                    ", the stored data is EOFs and eof_ts of " + var_name(nv)
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"neval","case","time","lat","lon"/)
dimSizes = (/neval  ,ncase ,ntime ,nlat ,nlon/)
dimUnlim = (/False  , False,False ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"case","integer",(/"case"/))
filevardef(ncdf,"time","integer",(/"time"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"eof",typeof(eof),(/"case","neval","lat","lon"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"case","neval","time"/))
filevardef(ncdf,"pcvar",typeof(pcvar),(/"case","neval"/))
filevardef(ncdf,"eval",typeof(eval),(/"case","neval"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"eof",eof)
filevarattdef(ncdf,"eof_ts",eof_ts)

ncdf->neval  = (/1,2,3,4/)  
ncdf->case   = (/1,2,3,4/)  ;the three case are (/"CTRL","NUDG"/)
ncdf->time   = ispan(1,360,1) 
ncdf->lat    = (/var&lat/)
ncdf->lon    = (/var&lon/)
ncdf->eof    = (/eof/)
ncdf->eof_ts = (/eof_ts/)
ncdf->pcvar  = (/pcvar/)
ncdf->eval   = (/eval/)
end

