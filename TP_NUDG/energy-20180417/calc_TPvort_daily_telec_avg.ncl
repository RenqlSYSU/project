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
fileout = path + "/energy-20180417/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-TPJJAvort_daily_telec_avg.nc"
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

f          = addfile(filein(0) + years(0) + ".daily." + var_name(0)+".nc","r")
time       = cd_calendar(f->time,0)
select_time= ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime      = dimsizes(select_time)
var  := f->$var_name(0)$(0,{lev},{flats:flatn},{flonl:flonr})
narea = dimsizes(area)
nyear = dimsizes(years)
ncase = dimsizes(case)
nflat = dimsizes(var&lat)
nflon = dimsizes(var&lon)
index = new((/narea,nyear,ntime/),float)
fvars = new((/nyear,ntime,nflat,nflon/),float)
cor   = new((/narea,nyear,nflat,nflon/),float)

do nf = 0,ncase-1,1
;==========================================================
;calculate the index firstly, year from 1979 to 2005
;============================================================
do nt = 0,nyear-1,1
    fu   = addfile(filein(nf) + years(0) + ".daily." + var_name(0) + ".nc","r")
    fv   = addfile(filein(nf) + years(0) + ".daily." + var_name(1) + ".nc","r")

    do na=0,narea-1,1
        u := fu->$var_name(0)$(select_time,{lev},{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
        v := fv->$var_name(1)$(select_time,{lev},{ilats(na):ilatn(na)},{ilonl(na):ilonr(na)})
        vort := uv2vr_cfd(u, v, v&lat, v&lon, 0) ;calculate the divergence
        index(na,nt,:) = wgt_areaave_Wrap(vort,1.0,1.0,0)  ; calculate the area average and get the index
    end do   ;the loop of the areas

    u  := fu->$var_name(0)$(select_time,{lev},{flats:flatn},{flonl:flonr})
    v  := fv->$var_name(1)$(select_time,{lev},{flats:flatn},{flonl:flonr})
    fvars(nt,:,:,:) = uv2vr_cfd(u, v, v&lat, v&lon, 0) ;calculate the divergence
    cor(:,nt,:,:)   = escorc_n(index(:,nt,:),fvars(nt,:,:,:),1,0)  ;get an array of (/narea,nyear,nlat,nlon/)
end do

corr   = dim_avg_n_Wrap(cor,1)
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

