;******************************************************************************
;filename:      calc_JJApreci_EOF.ncl
;output_file:   JJApreci_EOFn3.nc
;calculate the Asia JJA average precipitation (PRECC + PRECL) and 
;use the function of eofunc_n_Wrap, eofunc_ts_n_Wrap to 
;calculate the EOF pattern and time serises.
;                                            by Ql_Ren
;                                           2017/12/18
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-JJApreci_EOFn3.nc"
years    = ispan(1979,2005,1)   ;the calculate time
case     = (/"CTRL","NUDG24h","NUDG6h"/)
var_name = (/"PRECC","PRECL"/)

;-----area of east asia 
lats = 10
latn = 50
lonl = 110
lonr = 140
neval = 4

;------find the time index of the JJA and define the dimension of the varreas_telec_preci_JJA_corr.nc
f     = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars := f->$var_name(0)$(0,{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
ncase = dimsizes(case)
nvar  = dimsizes(var_name)
var   = new((/nyear,365,nlat,nlon/),float)
vard  = new((/nyear*365,nlat,nlon/),float)

do nf = 0,ncase-1,1 
;---------------------------------------------------------------
;read the data, then add the PRECC and PRECL
;--------------------------------------------------------------
opt = True
opt@jopt = 1
var = 0
do nt = 0,nyear-1,1
    do nv = 0,nvar-1,1 
        f  := addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
        var(nt,:,:,:) = var(nt,:,:,:) + f->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
    end do
end do  ;the loop of the years
var         = var*1000*3600*24
time       := cd_calendar(f->time,0)
select_time:= ind(time(:,1).ge.6.and.time(:,1).le.8)
avg_var     = dim_avg_n_Wrap(var(:,select_time,:,:),1) ;get the (nyear,nlat,nlon)
eof27       = eofunc_n_Wrap(avg_var,neval,opt,0) ;get the (neval,nlat,nlon)
eof27_ts    = eofunc_ts_n_Wrap(avg_var,eof27,False,0) ;get the (neval,nyear)
sig27       = eofunc_north(eof27@pcvar,nyear,True)
printVarSummary(eof27)
printVarSummary(eof27_ts)

vard = 0
do nv = 0, nvar-1, 1
    system_order = "ls " + filein(nf) + "*" + var_name(nv) + ".nc"
    files = systemfunc(system_order)
    f    := addfiles(files,"r")
    ListSetType(f,"cat")
    vard  = vard +  f[:]->$var_name(nv)$(:,{lats:latn},{lonl:lonr})
end do
vard        = vard*1000*3600*24
times       = f[:]->time
time       := cd_calendar(times,0)
select_time:= ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
eof27d      = eofunc_n_Wrap(vard(select_time,:,:),neval,opt,0) ;get the (neval,nlat,nlon)
eof27d_ts   = eofunc_ts_n_Wrap(vard(select_time,:,:),eof27d,False,0) ;get the (neval,nyear)
sig27d      = eofunc_north(eof27d@pcvar,ntime,True)
printVarSummary(eof27d)
printVarSummary(eof27d_ts)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "there are two types of EOF. The eof27 mean its sample number is 27 years" +\
                "and the eof27d mean its sample number is 27*92"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"neval","year","day","lat","lon"/)
dimSizes = (/neval  ,nyear ,ntime,nlat ,nlon /)
dimUnlim = (/False  ,False ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"day",typeof(times),(/"day"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"eof27"    ,typeof(eof27)    ,(/"neval","lat","lon"/))
filevardef(ncdf,"eof27_ts" ,typeof(eof27_ts) ,(/"neval","year"/))
filevardef(ncdf,"eof27d"   ,typeof(eof27d)   ,(/"neval","lat","lon"/))
filevardef(ncdf,"eof27d_ts",typeof(eof27d_ts),(/"neval","day"/))

;Define var Attribute
filevarattdef(ncdf,"day",times)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"eof27"    ,eof27)
filevarattdef(ncdf,"eof27_ts" ,eof27_ts)
filevarattdef(ncdf,"eof27d"   ,eof27d)
filevarattdef(ncdf,"eof27d_ts",eof27d_ts)

ncdf->neval  = (/1,2,3,4/)  
ncdf->year   = (/years/) 
ncdf->day    = (/times(select_time)/) 
ncdf->lat    = (/vars&lat/)
ncdf->lon    = (/vars&lon/)
ncdf->eof27     = (/eof27/)
ncdf->eof27_ts  = (/eof27_ts/)
ncdf->eof27d    = (/eof27d/)
ncdf->eof27d_ts = (/eof27d_ts/)
end do
end

