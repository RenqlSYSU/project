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
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-JJAvort_daily_EOF.nc"
years    = ispan(1979,2005,1)   ;the calculate time
case     = (/"CTRL","NUDG24h","NUDG6h"/)
var_name = (/"U","V"/)

;-----area of east asia 
lats = 25
latn = 40
lonl = 70
lonr = 130
lev  = 500 
neval = 8

;------find the time index of the JJA and define the dimension of the varreas_telec_preci_JJA_corr.nc
files1     = systemfunc("ls " + filein(0) + "*" + var_name(0) + ".nc")
f          = addfiles(files1,"r")
ListSetType(f,"cat")
time        = cd_calendar(f[:]->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
var  := f[:]->$var_name$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(var&lat)
nlon  = dimsizes(var&lon)
nyear = dimsizes(years)
ncase = dimsizes(case)

opt = True
opt@jopt = 1 ; uses correlation matrix to compute EOFs, 0 means use a covariance matrix 
do nf = 0,ncase-1,1 
;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
filu = systemfunc("ls " + filein(nf) + "*" + var_name(0) + ".nc")
fu   = addfiles(filu,"r")
ListSetType(fu,"cat")
u = fu[:]->$var_name(0)$(select_time,{lev},{lats:latn},{lonl:lonr})

filv  = systemfunc("ls " + filein(nf) + "*" + var_name(1) + ".nc")
fv    = addfiles(filv,"r")
ListSetType(fv,"cat")
v = fv[:]->$var_name(1)$(select_time,{lev},{lats:latn},{lonl:lonr})

vort  := uv2vr_cfd(u, v, v&lat, v&lon, 0) ;calculate the divergence
eof    = eofunc_n_Wrap(vort,neval,opt,0)         ;get the (neval,nlat,nlon)
eof_ts = eofunc_ts_n_Wrap(vort,eof,False,0) ;get the (neval,nyear)
sig    = eofunc_north(eof@pcvar,ntime,True)
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
dimNames = (/"neval","time","lat","lon"/)
dimSizes = (/neval  ,ntime ,nlat ,nlon/)
dimUnlim = (/False  ,False ,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"time" ,typeof(u&time),(/"time"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"eof",typeof(eof),(/"neval","lat","lon"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"neval","time"/))

;Define var Attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
filevarattdef(ncdf,"eof",eof)
filevarattdef(ncdf,"eof_ts",eof_ts)

ncdf->neval  = ispan(1,neval,1)  
ncdf->time   = (/u&time/) 
ncdf->lat    = (/var&lat/)
ncdf->lon    = (/var&lon/)
ncdf->eof    = (/eof/)
ncdf->eof_ts = (/eof_ts/)
end do
end

