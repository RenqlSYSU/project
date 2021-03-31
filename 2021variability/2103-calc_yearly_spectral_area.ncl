;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
var_name = (/"preci"/)
area = (/"EA","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 , 22  , 30  /)
latn = (/ 40 , 30  , 45  /)
lonl = (/105 ,107  ,115  /)
lonr = (/150 ,130  ,130  /)
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
nv = 0
na = 0
month_s = 5
month_e = 9

if(case_option.eq.0) then
    case    = (/"CTRL","NG15S","NG4550","NG45N","NG1520"/) ;
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/AMIP_NG15S/AMIP_CAM5_NG15S.cam.h1.YEAR.",\
                      "/model/AMIP_NG4550/AMIP_CAM5_NG4550.cam.h1.YEAR.",\
                      "/model/AMIP_NG45N/AMIP_CAM5_NG45N.cam.h1.YEAR.",\
                      "/model/AMIP_NG1520/AMIP_CAM5_NG1520.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
end if

do na = 0,0,1
fileout  = path + "/project/2021variability/data/" + case + "-" + month_s + month_e + "_spec_" + var_name(nv) + "-" + area(na) + ".nc"

f    = addfile(filein(0)+years(9)+".daily.PRECC.nc","r")
vars:= f->PRECC(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
if(month_s.gt.month_e) then 
select_time = ind(time(:,1).ge.month_s.or.time(:,1).le.month_e)
else
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
end if
ntime = dimsizes(select_time)
ncase = dimsizes(case)
nyear = dimsizes(years)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)

;---------------------------------------------------------------
;read the data, then calculate the area and season average
;--------------------------------------------------------------
; set spectral function arguments
  d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
  sm  = 1  ; smooth: should be at least 3 and odd; if small than 3, it would do no smoothing
  pct = 0;.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. If pct =0.0, no tapering will be done. If pct = 1.0, the whole series is affected

do nc = 0,2,1; ncase-1 ,1
print("handle with " + case(nc) + " " + area(na))
var  := new((/nyear,ntime,nlat,nlon/),float)
spcx := new((/nyear,4,toint(ntime/2),nlat,nlon/),float)

do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily.PRECC.nc","r")
    var(nt,:,:,:) = f->PRECC(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
    f   = addfile(filein(nc)+years(nt)+".daily.PRECL.nc","r")
    var(nt,:,:,:) = var(nt,:,:,:) + f->PRECL(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
var = var*1000*3600*24 ;convert m/s to mm/day
var = var - conform(var,dim_avg_n(var,0),(/1,2,3/))  ;remove annual cycle
printVarSummary(var)
var  = linmsg_n(var, -1, 3)
var  = linmsg_n(var, -1, 2)
printVarSummary(var)
if (any(ismissing(var))) then
    print("var has missing values")
end if

do nt = 0, nyear-1, 1
do ny = 0, nlat-1, 1
do nx = 0, nlon-1, 1
    ; calculate spectrum
    sdof = specx_anal(var(nt,:,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
    spcx(nt,:,:,ny,nx) = specx_ci(sdof, 0.05, 0.95) 
end do 
end do 
end do 
print(sdof)
printVarSummary(spcx)

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
fileatt@discription   = "yearly spectral powers, including input spectrum (indx=0), " + \
    "Markov Red Noise spectrum (indx=1), lower confidence bound for Markov (indx=2), " + \
    "upper confidence bound for Markov (indx=3)"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year","indx","lat","lon","frq"/);
dimSizes = (/nyear ,  4   ,nlat ,nlon ,nfrq /) ;
dimUnlim = (/False ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"frq" ,typeof(frq),(/"frq"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"spcx",typeof(spcx),(/"year","indx","frq","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"frq",frq)
filevarattdef(ncdf,"spcx",spcx)

ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->frq  = (/frq/)
ncdf->spcx = (/spcx/)
end do
end do
end
