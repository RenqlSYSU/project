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
lonl = 90 ;120
lonr = 300;260
;lev  = 200
neval = 4

opt = True
opt@jopt = 0 ;0 mean use the covariance matrix to compute EOFs (default), 1 mean use the correlation matrix
eof_kind = (/"cov","corr"/)
tname    = (/"_","_DJF_"/)

;----path of the file and output file-----------------------------
path     = "/home/ys17-19/renql/project/CONV_MON_WRK-2018"
gtopfile = path + "/data/USGS-gtopo30_1.9x2.5_remap_c050602.nc"
filename = path + "/data/" + case + "/pro/B20f19-" + (/"pure-aqua","surf-aqua","topo","realworld"/) + ".cam.h0." ; ; + var_name + ".1-100.nc" 
fileout  = path + "/mdata/" + var_name(nv) + lats + "S-" + latn + "N_" + lonl + "E-" + lonr + "E" + tname(0) + eof_kind(opt@jopt) + "EOF" + neval + ".nc"

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
var   = new((/ncase,nyear*nmonth,nlat,nlon/),float)
to_eof= new((/ncase,ntime,nlat,nlon/),float)
eof   = new((/ncase,neval,nlat,nlon/),float)
eof_ts= new((/ncase,neval,ntime/),float)
pcvar = new((/ncase,neval/),float)
eval  = new((/ncase,neval/),float)

f    = addfile(gtopfile,"r")
mlnd = f->LANDM_COSLAT({lats:latn},{lonl:lonr}) ;ocean (0), continent (1), transition (0-1)
mlnd_tmp = conform(var(0,:,:,:),mlnd,(/1,2/))
printVarSummary(mlnd_tmp)

;-----read the data-----------------------------
do nc = 0,ncase-1,1

if(ntime.eq.(nyear*nmonth)) then
;==========================================================================================
;apply EOF to all months(which have removal annual cycle) to find the evolution of ENSO 
;=======================dd=================================================================
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,:,:,:) = f->$var_name(nv)$(nd_s:nd_e,{lats:latn},{lonl:lonr})
    if(nc.ge.0) then
        var(nc,:,:,:) = mask(var(nc,:,:,:),mlnd_tmp.eq.1, False) ;remove land surface temperature
    end if
    var(nc,:,:,:) = var(nc,:,:,:) - 273.15  ;convert unit from K to C
    printVarSummary(var)
;remove the annual cycle
    to_eof(nc,:,:,:) = calcMonAnomTLL(var(nc,:,:,:),clmMonTLL(var(nc,:,:,:)))
else
;=============================================
;apply EOF to winter average (DJF) to find the dominant mode of mature ENSO
;==================================================
    f    = addfile(filename(nc) + var_name(nv) + ".1-100.nc","r")
    var(nc,:,:,:) = f->$var_name(nv)$((nd_s-1):(nd_e-1),{lats:latn},{lonl:lonr})
    if(nc.ge.2) then
        var(nc,:,:,:) = mask(var(nc,:,:,:),mlnd_tmp.eq.1, False) ;remove land surface temperature
    end if
    var(nc,:,:,:) = var(nc,:,:,:) - 273.15  ;convert unit from K to C
    printVarSummary(var)
;calculate the DJF average
    to_eof(nc,:,:,:) = (var(nc,0:359:12,:,:) + var(nc,1:359:12,:,:) + var(nc,2:359:12,:,:))/3
end if

    copy_VarMeta(var(0,0,:,:), to_eof(0,0,:,:))
    eof_temp    = eofunc_n_Wrap(to_eof(nc,:,:,:),neval,opt,0)         ;get the (neval,nlat,nlon)
    eof_ts_temp = eofunc_ts_n_Wrap(to_eof(nc,:,:,:),eof_temp,False,0) ;get the (neval,nyear)
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
ncdf->time   = ispan(1,ntime,1) 
ncdf->lat    = (/var&lat/)
ncdf->lon    = (/var&lon/)
ncdf->eof    = (/eof/)
ncdf->eof_ts = (/eof_ts/)
ncdf->pcvar  = (/pcvar/)
ncdf->eval   = (/eval/)
end

