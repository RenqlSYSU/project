;******************************************************************************
;filename:      calc_area_JJApreci_teleconnect.ncl
;output_file:   -Areas_JJApreci_telec.nc
;calculate the four areas JJA precipitation index and then 
;calculate the correlation between the index and other field precipitation 
;                                            by Ql_Ren
;                                           2017/12/02
;modification
;2018/01/24 calculate for three cases
;******************************************************************************
begin
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Areas_JJApreci_telec.nc"
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

f           = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
time        = cd_calendar(f->time,0)
var        := f->$var_name(0)$(:,{flats:flatn},{flonl:flonr})
select_time = ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime       = dimsizes(select_time)
nflat       = dimsizes(var&lat)
nflon       = dimsizes(var&lon)
nyear = dimsizes(years)
narea = dimsizes(area)
ncase = dimsizes(case)
index       = new((/narea,nyear/),float)
fvars       = new((/nyear,ntime,nflat,nflon/),float)
corr        = new((/narea,nflat,nflon/),float)

do nf = 0,ncase-1,1 
;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
do na=0,narea-1,1
    do t = 0,nyear-1,1
        index(na,t) = 0
        do v = 0, dimsizes(var_name)-1,1  ;the aim of this loop is to add the PRECC and PRECL 
            f  = addfile(filein(nf)+years(t)+".daily."+var_name(v)+".nc","r")
            vars := f->$var_name(v)$(select_time,{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
            vars := dim_avg_n_Wrap(vars,0)  ;calculate the JJA average
            index(na,t) = index(na,t) + wgt_areaave_Wrap(vars,1.0,1.0,0) ;calculate the area average and get the index
        end do ;the loop of precc and precl
    end do  ;the loop of the years
end do   ;the loop of the areas
printVarSummary(index)

;====================================================================
;read the field variable, then calculate the correlation between index and area
;======================================================================
do t = 0,nyear-1,1
    f_PRECC  = addfile(filein(nf)+years(t)+".daily."+var_name(0)+".nc","r")
    f_PRECL  = addfile(filein(nf)+years(t)+".daily."+var_name(1)+".nc","r")
    PRECC  = f_PRECC->$var_name(0)$(select_time,{flats:flatn},{flonl:flonr})
    PRECL  = f_PRECL->$var_name(1)$(select_time,{flats:flatn},{flonl:flonr})
    fvars(t,:,:,:) = PRECC + PRECL
end do
fvars1 = dim_avg_n_Wrap(fvars,1)  ;calculate the JJA average,the fvars is TLL
corr   = escorc_n(index,fvars1,1,0)  ;get an array of (/ncase,narea,nlat,nlon/) 
printVarSummary(fvars1)
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
                    "of four areas which are India,Tibet,Philippines,Maritime_continent"
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

