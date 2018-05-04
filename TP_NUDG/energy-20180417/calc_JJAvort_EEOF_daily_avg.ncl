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
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-JJAvort_daily_EEOF_avg.nc"
years    = ispan(1979,2005,1)   ;the calculate time
case     = (/"CTRL","NUDG24h","NUDG6h"/)
var_name = (/"U","V"/)

;-----area of east asia 
lats = 25
latn = 40
lonl = 70
lonr = 130
lev  = 500 
neval = 4
delay = 8

;------find the time index of the JJA and define the dimension of the varreas_telec_preci_JJA_corr.nc
f          = addfile(filein(0) + years(0) + ".daily." + var_name(0)+".nc","r")
time       = cd_calendar(f->time,0)
select_time= ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime      = dimsizes(select_time)
var  := f->$var_name$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(var&lat)
nlon  = dimsizes(var&lon)
nyear = dimsizes(years)
ncase = dimsizes(case)
u     = new((/nyear,ntime,nlat,nlon/),float)
v     = new((/nyear,ntime,nlat,nlon/),float)
vort  = new((/nyear,ntime,nlat,nlon/),float)
vort_temp   = new((/nyear,delay,(ntime-delay),nlat,nlon/),float)
eof_temp    = new((/nyear,neval,delay,nlat,nlon/),float)
eof_ts_temp = new((/nyear,neval,ntime-delay/),float)

opt = True
opt@jopt = 1 ; uses correlation matrix to compute EOFs, 0 means use a covariance matrix 
do nf = 0,ncase-1,1 
;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
do nt = 0,nyear-1,1
    fu   = addfile(filein(nf) + years(0) + ".daily." + var_name(0) + ".nc","r")
    fv   = addfile(filein(nf) + years(0) + ".daily." + var_name(1) + ".nc","r")
    u(nt,:,:,:) = fu->$var_name(0)$(select_time,{lev},{lats:latn},{lonl:lonr})
    v(nt,:,:,:) = fv->$var_name(1)$(select_time,{lev},{lats:latn},{lonl:lonr})
    vort(nt,:,:,:)  = uv2vr_cfd(u(nt,:,:,:), v(nt,:,:,:), v&lat, v&lon, 0) ;calculate the divergence
end do

do ndelay = 0,delay-1,1
    time_end = ntime - delay + ndelay-1
    vort_temp(:,ndelay,:,:,:) = vort(:,ndelay:time_end,:,:)
end do

do nt = 0,nyear-1,1
    eof_temp(nt,:,:,:,:) = eofunc_n_Wrap(vort_temp(nt,:,:,:,:),neval,opt,1)         ;get the (neval,nlat,nlon)
    eof_ts_temp(nt,:,:)  = eofunc_ts_n_Wrap(vort_temp(nt,:,:,:,:),eof_temp(nt,:,:,:,:),False,1) ;get the (neval,nyear)
    sig    := eofunc_north(eof_temp@pcvar,ntime,True)
end do
eof    = dim_avg_n_Wrap(eof_temp,0)   ;get (delay,neval,nlat,nlon)
eof_ts = dim_avg_n_Wrap(eof_ts_temp,0)
printVarSummary(eof)
printVarSummary(eof_ts)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the stored data is EOFs and eof_ts of JJA vort."
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"neval","delay","time"     ,"lat","lon"/)
dimSizes = (/neval  ,delay  ,ntime-delay,nlat ,nlon/)
dimUnlim = (/False  ,False  ,False      ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"delay","integer",(/"delay"/))
filevardef(ncdf,"time" ,typeof(u&time),(/"time"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"eof",typeof(eof),(/"neval","delay","lat","lon"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"neval","time"/))

;Define var Attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"eof",eof)
filevarattdef(ncdf,"eof_ts",eof_ts)

ncdf->neval  = ispan(1,neval,1)  
ncdf->delay  = ispan(1,delay,1)  
ncdf->time   = (/u&time(0:ntime-delay-1)/) 
ncdf->lat    = (/var&lat/)
ncdf->lon    = (/var&lon/)
ncdf->eof    = (/eof/)
ncdf->eof_ts = (/eof_ts/)
end do
end

