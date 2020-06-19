;***************************************************************************
;use the f_Qd, f_Qd_t, f_Qeddy and A to calc dzdt by SOR (Successive Over Relaxation)
;1. read f_Qd, f_Qd_t, f_Qeddy and A
;2. define coefficient(coe) of linear equations and relaxing factor(rf)
;3. calc dzdt by circulation 
;4. store the results
;
;The ncl functions have been used are: 
;   dim_avg_n_Wrap, center_finite_diff_n, conform
;   dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/10/25
;******************************************************************************
begin
path   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
data_type = (/"double","float"/)
ntype = 0
ntime = 37 
freq = 4 

case   = "ERA-Interim"
;case   = "NCEP1"
if(freq.eq.0) then 
filein  = path + case + "_DJF_Qeddy_A_freq-year.nc"
filedat = path + case + "_DJF_dzdt_freq.dat"
fileout = path + case + "_DJF_dzdt_freq-year.nc"
var_name= (/"f_Qeddh","f_Qeddl","Ah","Al"/)  
end if
if(freq.eq.1) then 
filein  = path + case + "_DJF_Qd_Qeddy_A-year.nc"
filedat = path + case + "_DJF_dzdt.dat"
fileout = path + case + "_DJF_dzdt-year.nc"
var_name= (/"f_Qd","f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
end if
if(freq.eq.2) then 
filein  = path + case + "_DJF_forc_in-year.nc"
filedat = path + case + "_DJF_dzdt_advect.dat"
fileout = path + case + "_DJF_dzdt_advect-year.nc"
var_name= (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
end if
if(freq.eq.3) then 
filein  = path + case + "_DJF_forc_in2-year.nc"
filedat = path + case + "_DJF_dzdt_advect_xy.dat"
fileout = path + case + "_DJF_dzdt_in3-year.nc"
var_name= (/"b3x","b3y","b1x","b1y"/)  ;this is the variable used to draw, can used as center string
end if
if(freq.eq.4) then 
filein  = path + case + "_DJF_forc_re-year.nc"
filedat = path + case + "_DJF_dzdt_re.dat"
fileout = path + case + "_DJF_dzdt_re-year.nc"
var_name= (/"turn","b1z","b4"/)  ;this is the variable used to draw, can used as center string
end if

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats = 10 ;-20 ;-90 ;
latn = 90  ; 90 ;
lonl =  0  ;30  ;
lonr = 360 ;150 ;
lev  = (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,175,150,125, 100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein,"r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)
if(ntype.eq.0) then 
dzdt  = new((/nvar,ntime,nlev,nlat,nlon/),double)
else
dzdt  = new((/nvar,ntime,nlev,nlat,nlon/),float)
end if

itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
coslat = conform(dzdt,cos(lat),ilat)

dlev = lev
dlev(1:(nlev-1)) = (- lev(0:(nlev-2)) + lev(1:(nlev-1)))*100 
dlev(0) = dlev(1) 
dlev5  = conform(dzdt,dlev,ilev)

opt    = 0     ;used by center_finite_diff_n, no meanging 
loncyc = False;True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
setfileoption("bin","ReadByteOrder","LittleEndian")
irec = 0
do nv = 0, nvar-1, 1
do nt = 0, ntime-1, 1
do nz = 0, nlev-1, 1
    dzdt(nv,nt,nz,:,:) = fbindirread(filedat,irec,(/nlat,nlon/),data_type(ntype))
    irec = irec + 1
end do
end do
end do
print(dzdt(0:(nvar-1),0,15,29,119))
print(dzdt(0,0,15,29,129:149:4))
print(dzdt(0,0,9 ,29,129:149:4))

;==========================================================
;define coefficient and calc dzdt
;============================================================
if(ntype.eq.0) then 
f0   = todouble(conform(dzdt,2*(2*pi/24.0/3600.0)*sin(lat),ilat))
else
f0   = tofloat(conform(dzdt,2*(2*pi/24.0/3600.0)*sin(lat),ilat))
end if
f0@_FillValue = dzdt@_FillValue
printVarSummary(f0)

do ny = 0, nlat-1
if (abs(lat(ny)*180.0/pi ).lt. 10. ) then
f0(:,:,:,ny,:)= f0@_FillValue
end if
end do
printVarSummary(dzdt)
dudt = - center_finite_diff_n(dzdt,lat,cyclic,opt,ilat)/f0/a
dvdt =   center_finite_diff_n(dzdt,lon,loncyc,opt,ilon)/f0/a/coslat
dTdt = - center_finite_diff_n(dzdt,lev*100,cyclic,opt,ilev)*conform(dzdt,lev*100,ilev)/R
dbdt = - center_finite_diff_n(dTdt,lat,cyclic,opt,ilat)/a

;nvar@long_name = "dzdt induced by f_Qd, f_Qd_t, f_Qeddy and A" 
dzdt@long_name = "Climatology daily geopotential tendency"
dudt@long_name = "Climatology daily zonal wind tendency"
dvdt@long_name = "Climatology daily meridional wind tendency"
dTdt@long_name = "Climatology daily temperature tendency"
dbdt@long_name = "Climatology daily baroclinicity(-dTdy) tendency"
dzdt@units     = "m^2/s^3" 
dudt@units     = "m/s^2" 
dvdt@units     = "m/s^2" 
dTdt@units     = "K/s" 
dbdt@units     = "K/(m*s)" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "DJF"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"var","year","lev","lat","lon"/);
dimSizes = (/nvar ,ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"nvar","integer",(/"var"/))
filevardef(ncdf,"lev" ,typeof(vars&lev),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"dzdt",typeof(dzdt),(/"var","year","lev","lat","lon"/))
filevardef(ncdf,"dudt",typeof(dudt),(/"var","year","lev","lat","lon"/))
filevardef(ncdf,"dvdt",typeof(dvdt),(/"var","year","lev","lat","lon"/))
filevardef(ncdf,"dTdt",typeof(dTdt),(/"var","year","lev","lat","lon"/))
filevardef(ncdf,"dbdt",typeof(dbdt),(/"var","year","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"lev",vars&lev)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"dzdt",dzdt)
filevarattdef(ncdf,"dudt",dudt)
filevarattdef(ncdf,"dvdt",dvdt)
filevarattdef(ncdf,"dTdt",dTdt)
filevarattdef(ncdf,"dbdt",dbdt)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->dzdt = (/dzdt/)
ncdf->dudt = (/dudt/)
ncdf->dvdt = (/dvdt/)
ncdf->dTdt = (/dTdt/)
ncdf->dbdt = (/dbdt/)
end


