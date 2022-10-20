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
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-TPJJAvort_daily_telec.nc"
var_name = (/"U","V"/)
case     = (/"CTRL","NG24h","NG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time

;area of index, India,Tibet,Philippines,Maritime_continent
area  = (/"CTP"/)
ilats = (/ 30/)
ilatn = (/ 35/)
ilonl = (/ 80/)
ilonr = (/100/)
lev = (/500/)

;area of field variables,larger then the plan
flats = (/-20/)
flatn = (/60 /)
flonl = (/25 /)
flonr = (/155/)

files1     = systemfunc("ls " + filein(0) + "*" + var_name(0) + ".nc")
f          = addfiles(files1,"r")
ListSetType(f,"cat")
time       = cd_calendar(f[:]->time,0)
select_time= ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime      = dimsizes(select_time)
var  := f[:]->$var_name(0)$(0,{lev},{flats:flatn},{flonl:flonr})
narea = dimsizes(area)
ncase = dimsizes(case)
nflat = dimsizes(var&lat)
nflon = dimsizes(var&lon)
index = new((/narea,ntime/),float)

do nf = 0,ncase-1,1
;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
    filu = systemfunc("ls " + filein(nf) + "*" + var_name(0) + ".nc")
    fu   = addfiles(filu,"r")
    ListSetType(fu,"cat")
    filv  = systemfunc("ls " + filein(nf) + "*" + var_name(1) + ".nc")
    fv    = addfiles(filv,"r")
    ListSetType(fv,"cat")

    do na=0,narea-1,1
        u := fu[:]->$var_name(0)$(select_time,{lev},{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
        v := fv[:]->$var_name(1)$(select_time,{lev},{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
        vort := uv2vr_cfd(u, v, v&lat, v&lon, 0) ;calculate the divergence
        index(na,:) = wgt_areaave_Wrap(vort,1.0,1.0,0)  ; calculate the area average and get the index
    end do   ;the loop of the areas

    u  := fu[:]->$var_name(0)$(select_time,{lev},{flats:flatn},{flonl:flonr})
    v  := fv[:]->$var_name(1)$(select_time,{lev},{flats:flatn},{flonl:flonr})
    fvars = uv2vr_cfd(u, v, v&lat, v&lon, 0) ;calculate the divergence

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
fileAtt@discription   = "the correlation between Asia JJA vorticity and area average JJA vort index ,area is TP"
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
corr@long_name = "the teleconnection of Aisa with vort index in JJA"
corr@units     = ""
filevarattdef(ncdf,"corr",corr)

ncdf->area= (/1/) ;the four areas are (/"India","Tibet","Philippines","Maritime_continent"/)
ncdf->lat = (/var&lat/)
ncdf->lon = (/var&lon/)
ncdf->corr= (/corr/)
end do
end

