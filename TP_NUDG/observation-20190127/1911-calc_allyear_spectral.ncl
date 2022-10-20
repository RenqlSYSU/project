;***************************************************************************
;filename:      spec.ncl
;read the 27 years U,V,T at TP areas in three cases (CTRL, NUDG6h, NUDG24h)
;then do the spectral analysis and plot the figure
;                                           by Ql_Ren
;                                           2019/01/20
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
var_name = (/"uwnd","vwnd","air","vort"/)
area = (/"EA","TP","SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 25 ,25  , 22  , 30  /)
latn = (/ 40 ,40  , 30  , 45  /)
lonl = (/ 0  ,70  ,107  ,115  /)
lonr = (/360 ,105 ,130  ,130  /)
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
nv = 1
na = 0
area_option = True

if(na.eq.0) then
lev  = (/850,500,400,300,200/) ;850,500,
if(area_option) then 
fileout  = path2 + "mdata/ERA-Interim_allyear_spec_" + var_name(nv) + "-" + area(na) + "-area.nc"
else
fileout  = path2 + "mdata/ERA-Interim_allyear_spec_" + var_name(nv) + "-" + area(na) + ".nc"
end if
else
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
end if
years    = ispan(1979,2005,1)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars:= f->$var_name(0)$(0,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
nlev  = dimsizes(lev)
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

print("handle with ERA-Interim " + area(na))
if(nv.eq.3) then 
    inv = 0
    files = systemfunc("ls " + path1 + var_name(inv) + "/" + var_name(inv) + ".interim.*.nc")
    print(files(0:nyear-1))
    f := addfiles(files(0:nyear-1),"r")
    ListSetType(f,"cat")
    u = f[:]->$var_name(inv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    inv = 1
    files = systemfunc("ls " + path1 + var_name(inv) + "/" + var_name(inv) + ".interim.*.nc")
    print(files(0:nyear-1))
    f := addfiles(files(0:nyear-1),"r")
    ListSetType(f,"cat")
    v = f[:]->$var_name(inv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    var = uv2vr_cfd(u,v,vars&lat,vars&lon,2)
    copy_VarMeta(u,var)
    delete([/u,v/])
else
    files = systemfunc("ls " + path1 + var_name(nv) + "/" + var_name(nv) + ".interim.*.nc")
    print(files(0:nyear-1))
    f := addfiles(files(0:nyear-1),"r")
    ListSetType(f,"cat")
    var = f[:]->$var_name(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end if
printVarSummary(var)
if (any(ismissing(var))) then
    print("var1 has missing values")
end if

ntime = dimsizes(var&time)
spcx := new((/toint(ntime/2),nlev,nlat,nlon/),float)
do nz = 0, nlev-1, 1
do ny = 0, nlat-1, 1
do nx = 0, nlon-1, 1
; calculate spectrum
    sdof = specx_anal(var(:,nz,ny,nx),d,sm,pct) ;returns the degrees of freedom as a scalar.
    spcx(:,nz,ny,nx) = sdof@spcx 
end do 
end do 
end do 
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
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
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
end
