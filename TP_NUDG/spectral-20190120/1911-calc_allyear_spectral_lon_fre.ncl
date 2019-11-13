;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
var_name = (/"U","V","T"/)
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case    = (/"CTRL","NUDG6h","NUDG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
end if

if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case1   = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    case    = (/"CTRL","noTP"/)
    filein  = path + "/model/" + case1 + "/pro/" + case1 + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
end if

if(case_option.eq.2) then
    case    = (/"F2000","NG6h"/)
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_NUDG/F2000_CAM5.cam.h1.YEAR.00"/)
    years    = ispan(1,30,1)   ;the calculate time
end if
area = (/"EA","TP","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 ,25  , 22  , 30  /)
latn = (/ 40 ,40  , 30  , 45  /)
lonl = (/ 30 ,70  ,107  ,115  /)
lonr = (/150 ,105 ,130  ,130  /)
nv = 1
na = 1
if(na.eq.0) then
lev  = (/850,500,200/) ;850,500,
fileout  = path + "/project/TP_NUDG/spectral-20190120/mdata/" + case + "-allyear_spec_" + var_name(nv) + "-" + area(na) + ".nc"
else
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
;lev  = (/850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
fileout  = path + "/project/TP_NUDG/spectral-20190120/mdata/" + case + "-allyear_spec_" + var_name(nv) + "-" + area(na) + ".nc"
end if

f    = addfile(filein(0)+years(9)+".daily."+var_name(nv)+".nc","r")
vars:= f->$var_name(nv)$(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
altime = ntime*nyear

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0,1,1; ncase-1 ,1
print("handle with " + case(nc) + " " + area(na))
var  := new((/nyear,ntime,nlev,nlat,nlon/),float)
spcx := new((/toint(altime/2),nlev,nlat,nlon/),float)

if((case_option.eq.1).and.(nc.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
if((case_option.eq.1).and.(nc.eq.0))
    years   = ispan(21,30,1)   ;the calculate time
end if
files = systemfunc("ls "+filein(nc)+"*.daily."+var_name(nv)+".nc" )
do nt = 0, nyear-1 ,1
    if(case_option.eq.2) then 
        f   = addfile(files(nt),"r")
    else
        f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    end if
    var(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
;var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle
printVarSummary(var)
var  = linmsg_n(var, -1, 4)
var  = linmsg_n(var, -1, 3)
var  = linmsg_n(var, -1, 2)
var1 = reshape(var,(/altime,nlev,nlat,nlon/))
printVarSummary(var)
if (any(ismissing(var1))) then
    print("var1 has missing values")
end if

do nz = 0, nlev-1, 1
do ny = 0, nlat-1, 1
do nx = 0, nlon-1, 1
; calculate spectrum
    sdof = specx_anal(var1(:,nz,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
    spcx(:,nz,ny,nx) = sdof@spcx 
end do 
end do 
end do 
print(sdof)
printVarSummary(spcx)

if(na.eq.0) then 
    spcx := dim_avg_n(spcx, 2) ;get var(nvar,ncase,nyear,365,nlev)
else
    spcx := wgt_areaave(spcx, 1.0, 1.0, 0) ;get var(nvar,ncase,nyear,365,nlev)
end if
frq  = sdof@frq          ; convenience
nfrq = dimsizes(frq)
printVarSummary(spcx)
print(frq(0:50:10))

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "The frequency-longitude or level-frequency distribution of spectral powers"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"lev","lon","frq"/);
dimSizes = (/nlev ,nlon ,nfrq /) ;
dimUnlim = (/False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"frq" ,typeof(frq),(/"frq"/))
if(na.eq.0) then 
filevardef(ncdf,"spcx",typeof(spcx),(/"frq","lev","lon"/))
else
filevardef(ncdf,"spcx",typeof(spcx),(/"frq","lev"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"frq",frq)
filevarattdef(ncdf,"spcx",spcx)

ncdf->lev  = (/lev/)
ncdf->lon  = (/vars&lon/)
ncdf->frq  = (/frq/)
ncdf->spcx = (/spcx/)
end do
end
