;***************************************************************************
;filename:      calc_qtrans_ano_850hpa.ncl
;output_file:   CTRL-JJA_ano_850hpa_qu.nc 
;read the 850hPa 27 years JJA daily variable(TLL) of one case at East AsiaI(15S£­55N£¬30-150E) 
;and then calculate the climate average(LL) and the anomaly(TLL)
;                                            by Ql_Ren
;                                           2018/04/21
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-JJA_ano_850hpa_"
var_name  = (/"U","Q"/)
calc_name = (/"ave_qu","aveq_aveu","aveq_u","q_aveu","qu"/)
;var_name  = (/"V","Q"/)
;calc_name = (/"ave_qv","aveq_avev","aveq_v","q_avev","qv"/)
years    = ispan(1979,2005,1)   ;the calculate time
g  = 9.8 ;m/(s*s)

;area: east Asia
lats = -15 ;-20
latn = 55  ;60
lonl = 30  ;20
lonr = 150 ;220
lev  = (/850/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var_u      = new((/nyear,ntime,nlat,nlon/),float)
var_q      = new((/nyear,ntime,nlat,nlon/),float)

do nf = 0 , ncase-1,1
;---------------------------------------------------------------
;read the data
;--------------------------------------------------------------
nv = 0
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_u(nt,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years
printVarSummary(var_u)

nv = 1
do nt = 0,nyear-1,1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var_q(nt,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do   ;the loop of 27 years
var_q = var_q*1000  ;convert unit from kg/kg to g/kg
printVarSummary(var_q)

var = var_u*var_q/g  ;get var(nyear,ntime,nlat,nlon)
ave_qu = dim_avg_n_Wrap(var,(/0,1/)) ;get var(nlat,nlon)
ave_qu@long_name = "the average of 27 years and 92 days of 850hpa qu"
ave_qu@units     = "g/(s*pa*m)"

aveq = dim_avg_n_Wrap(var_q,0) ;get (ntime,nalt,nlon)
aveu = dim_avg_n_Wrap(var_u,0) ;get (ntime,nlat,nlon)
q    = dim_rmvmean_n_Wrap(var_q,0) ;get (nyear,ntime,nlat,nlon)
u    = dim_rmvmean_n_Wrap(var_u,0) ;get (nyear,ntime,nlat,nlon)

aveq_aveu = dim_avg_n_Wrap(aveq*aveu/g,0) ;get (nlat,nlon) 
aveq_u    = dim_avg_n_Wrap(conform(u,aveq,(/1,2,3/))*u/g,(/0,1/)) ;get (nlat,nlon)
q_aveu    = dim_avg_n_Wrap(conform(q,aveu,(/1,2,3/))*q/g,(/0,1/)) ;get (nlat,nlon)
qu        = dim_avg_n_Wrap(q*u,(/0,1/)) ;get (nlat,nlon)

;====================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+calc_name(4)+".nc")
ncdf = addfile(fileout(nf)+calc_name(4)+".nc","c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the climatology and JJA seasonal 850 hPa moisture transport analysis "
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"lat","lon"/)
dimSizes = (/nlat ,nlon/)
dimUnlim = (/False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"ave_qu"   ,typeof(ave_qu   ) ,(/"lat","lon"/))
filevardef(ncdf,"aveq_aveu",typeof(aveq_aveu) ,(/"lat","lon"/))
filevardef(ncdf,"aveq_u"   ,typeof(aveq_u   ) ,(/"lat","lon"/))
filevardef(ncdf,"q_aveu"   ,typeof(q_aveu   ) ,(/"lat","lon"/))
filevardef(ncdf,"qu"       ,typeof(qu       ) ,(/"lat","lon"/))

;Define the attribute of the variables
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"ave_qu",ave_qu)

ncdf->lat = (/vars&lat/)
ncdf->lon = (/vars&lon/)
ncdf->ave_qu    = (/ave_qu   /)
ncdf->aveq_aveu = (/aveq_aveu/)
ncdf->aveq_u    = (/aveq_u   /)
ncdf->q_aveu    = (/q_aveu   /)
ncdf->qu        = (/qu       /)
end do
end

