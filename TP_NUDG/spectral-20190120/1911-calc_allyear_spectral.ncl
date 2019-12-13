;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
var_name = (/"U","V","T","vort"/)
area = (/"EA","TP","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 ,25  , 22  , 30  /)
latn = (/ 40 ,40  , 30  , 45  /)
lonl = (/ 0  ,70  ,107  ,115  /)
lonr = (/360 ,105 ,130  ,130  /)
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
nv = 3
na = 0
month_s = 4
month_e = 9
area_option = False

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

do na = 1,3,1
if(na.eq.0) then
lev  = (/850,500,400,300,200/) ;850,500,
if(area_option) then 
fileout  = path + "/project/TP_NUDG/spectral-20190120/mdata/" + case + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + "-area.nc"
else
fileout  = path + "/project/TP_NUDG/spectral-20190120/mdata/" + case + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"
end if
else
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
;lev  = (/850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
fileout  = path + "/project/TP_NUDG/spectral-20190120/mdata/" + case + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"
end if

f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
if(ntime.eq.365) then 
altime = ntime*nyear
else
altime = ntime
end if

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0,2,1; ncase-1 ,1
print("handle with " + case(nc) + " " + area(na))
var  := new((/nyear,ntime,nlev,nlat,nlon/),float)
spcx := new((/toint(altime/2),nlev,nlat,nlon/),float)

if((case_option.eq.1).and.(nc.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
if((case_option.eq.1).and.(nc.eq.0))
    years   = ispan(21,30,1)   ;the calculate time
end if

if(nv.eq.3) then 
u = var
v = var
inv = 0
files = systemfunc("ls "+filein(nc)+"*.daily."+var_name(inv)+".nc" )
do nt = 0, nyear-1 ,1
    if(case_option.eq.2) then 
        f   = addfile(files(nt),"r")
    else
        f   = addfile(filein(nc)+years(nt)+".daily."+var_name(inv)+".nc","r")
    end if
    u(nt,:,:,:,:) = f->$var_name(inv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
inv = 1
files = systemfunc("ls "+filein(nc)+"*.daily."+var_name(inv)+".nc" )
do nt = 0, nyear-1 ,1
    if(case_option.eq.2) then 
        f   = addfile(files(nt),"r")
    else
        f   = addfile(filein(nc)+years(nt)+".daily."+var_name(inv)+".nc","r")
    end if
    v(nt,:,:,:,:) = f->$var_name(inv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
v  = linmsg_n(v, -1, 4)
v  = linmsg_n(v, -1, 3)
v  = linmsg_n(v, -1, 2)
u  = linmsg_n(u, -1, 4)
u  = linmsg_n(u, -1, 3)
u  = linmsg_n(u, -1, 2)
var = uv2vr_cfd(u,v,vars&lat,vars&lon,2)
delete([/u,v/])
else
files = systemfunc("ls "+filein(nc)+"*.daily."+var_name(nv)+".nc" )
do nt = 0, nyear-1 ,1
    if(case_option.eq.2) then 
        f   = addfile(files(nt),"r")
    else
        f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    end if
    var(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
end if

;var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle
printVarSummary(var)
var  = linmsg_n(var, -1, 4)
var  = linmsg_n(var, -1, 3)
var  = linmsg_n(var, -1, 2)
printVarSummary(var)
if (any(ismissing(var))) then
    print("var has missing values")
end if
;var1 = dim_standardize_n_Wrap(var1,1,0)

if(ntime.eq.365) then 
    var1 := reshape(var,(/altime,nlev,nlat,nlon/))
    do nz = 0, nlev-1, 1
    do ny = 0, nlat-1, 1
    do nx = 0, nlon-1, 1
        ; calculate spectrum
        sdof = specx_anal(var1(:,nz,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
        spcx(:,nz,ny,nx) = sdof@spcx 
    end do 
    end do 
    end do 
else
    spcx = 0.0
    do nt = 0, nyear-1, 1
    do nz = 0, nlev-1, 1
    do ny = 0, nlat-1, 1
    do nx = 0, nlon-1, 1
        ; calculate spectrum
        sdof = specx_anal(var(nt,:,nz,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
        spcx(:,nz,ny,nx) = spcx(:,nz,ny,nx) + sdof@spcx 
    end do 
    end do 
    end do 
    end do 
    spcx = spcx/tofloat(nyear)
end if
print(sdof)
printVarSummary(spcx)

if(na.eq.0) then 
    if(.not.(area_option)) then
    spcx := dim_avg_n(spcx, 2) ;get var(nvar,ncase,nyear,365,nlev)
    end if
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
dimNames = (/"lev","lat","lon","frq"/);
dimSizes = (/nlev ,nlat ,nlon ,nfrq /) ;
dimUnlim = (/False,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"frq" ,typeof(frq),(/"frq"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
if(na.eq.0) then 
if(area_option) then
filevardef(ncdf,"spcx",typeof(spcx),(/"frq","lev","lat","lon"/))
else
filevardef(ncdf,"spcx",typeof(spcx),(/"frq","lev","lon"/))
end if
else
filevardef(ncdf,"spcx",typeof(spcx),(/"frq","lev"/))
end if

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"frq",frq)
filevarattdef(ncdf,"spcx",spcx)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->frq  = (/frq/)
ncdf->spcx = (/spcx/)
end do
end do
end
