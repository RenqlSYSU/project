;******************************************************************************
;filename:      calc_area_JJApreci_daliy_telec.ncl
;output_file:   -Areas_JJApreci_daily_telec.nc
;calculate the four areas JJA precipitation index and then 
;calculate the correlation between the index and other field precipitation
;the sample number is 27years*92
;                                            by Ql_Ren
;                                           2017/01/24
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Areas_JJApreci_daily_telec.nc"
var_name = (/"PRECC","PRECL"/)
case     = (/"CTRL","NUDG24h","NUDG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"India","Tibet","Philippines","Maritime_continent"/)
ilats = (/20,25,5,-10/)
ilatn = (/35,35,20,7.5/)
ilonl = (/60,90,120,85/)
ilonr = (/80,103,140,120/)

;area of field variables,larger then the plan
flats = -20 ;-15
flatn = 60 ;55
flonl = 20 ;30
flonr = 160 ;150

files1     = systemfunc("ls " + filein(0) + "*" + var_name(0) + ".nc")
f_precc    = addfiles(files1,"r")
ListSetType(f_precc,"cat")
time        = cd_calendar(f_precc[:]->time,0)
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
var  := f_precc[:]->$var_name(0)$(0,{flats:flatn},{flonl:flonr})
narea = dimsizes(area)
ncase = dimsizes(case)
nflat = dimsizes(var&lat)
nflon = dimsizes(var&lon)
index       = new((/narea,ntime/),float)

do nf = 0,ncase-1,1 
;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
files1     = systemfunc("ls " + filein(nf) + "*" + var_name(0) + ".nc")
f_precc    = addfiles(files1,"r")
ListSetType(f_precc,"cat")
files2     = systemfunc("ls " + filein(nf) + "*" + var_name(1) + ".nc")
f_precl    = addfiles(files2,"r")
ListSetType(f_precl,"cat")

do na=0,narea-1,1
    precc := f_precc[:]->$var_name(0)$(select_time,{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
    precl := f_precl[:]->$var_name(1)$(select_time,{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
    index(na,:) = wgt_areaave_Wrap(precl,1.0,1.0,0) + wgt_areaave_Wrap(precc,1.0,1.0,0) ;calculate the area average and get the index
end do   ;the loop of the areas

precc  := f_precc[:]->$var_name(0)$(select_time,{flats:flatn},{flonl:flonr})
precl  := f_precl[:]->$var_name(1)$(select_time,{flats:flatn},{flonl:flonr})
fvars   = precc + precl              ;get an array of (/ntime,nlat,nlon/)
corr    = escorc_n(index,fvars,1,0)  ;get an array of (/narea,nlat,nlon/) 
printVarSummary(ntime)
printVarSummary(index)
printVarSummary(fvars)
printVarSummary(corr)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the correlation between Asia JJA rainfall and area average JJA rainfall index " + \
                    "of four areas which are India,Tibet,Philippines,Maritime_continent,and the sample number is 27*92"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"area","lat" ,"lon"/)
dimSizes = (/narea ,nflat ,nflon/)
dimUnlim = (/False ,False ,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"area","integer",(/"area"/))
filevardef(ncdf,"lat",typeof(var&lat),getvardims(var&lat))
filevardef(ncdf,"lon",typeof(var&lon),getvardims(var&lon))
filevardef(ncdf,"corr",typeof(corr),(/"area","lat","lon"/))

;Define var Attribute
filevarattdef(ncdf,"lat",var&lat)
filevarattdef(ncdf,"lon",var&lon)
corr@long_name = "the teleconnection of Aisa with preci (PRECC+PRECL) index in JJA"
corr@units     = ""
filevarattdef(ncdf,"corr",corr)

ncdf->area= (/1,2,3,4/) ;the four areas are (/"India","Tibet","Philippines","Maritime_continent"/)
ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr= (/corr/)
end do
end

