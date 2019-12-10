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
path    = "/home/ys17-19/renql"
obs_option = False
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
data_type = (/"double","float"/)
ntype = 0
ntime = 4
nc = 1

if(case_option.eq.0) then
    case   = (/"CTRL","NUDG6h","NUDG24h"/)
else
    case   = (/"TP_CTRL","TP_CR"/)
end if

if(.not.(obs_option)) then 
filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_Qd_Qeddy_A_lon.nc"
filedat = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_dzdt_month.dat"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_dzdt_month.nc"
else
if(ntime.eq.1) then 
filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_Qd_Qeddy_A_1204.nc"
filedat = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_dzdt.dat"
fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_dzdt.nc" 
else
filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_month_Qd_Qeddy_A2.nc"
filedat = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_dzdt_month.dat"
fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_dzdt_month.nc" 
end if
end if
;var_name= (/"f_Qd","f_Qd_t","f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
var_name= (/"f_Qd","f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
;var_name = (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats = 10 ;-20 ;-90 ;
latn = 90  ; 90 ;
lonl = 0  ; 0  ;
lonr = 360 ;360 ;
;lev  = (/1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100/)
;lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
;lev  = (/1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150, 100/) 
lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/)
;lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)
month_s = 1
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(nc),"r")
vars:= f->Qd(0,{lev},{lats:latn},{lonl:lonr})
;time = cd_calendar(f->time,0)
;select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
;ntime = dimsizes(select_time)
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
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 0, nc, 1 
print("handle with " + case(nf))
;==========================================================
;read the data
;============================================================
setfileoption("bin","ReadByteOrder","LittleEndian")
irec = 0
do nv = 0, nvar-1, 1
do nt = 0, ntime-1, 1
do nz = 0, nlev-1, 1
    dzdt(nv,nt,nz,:,:) = fbindirread(filedat(nf),irec,(/nlat,nlon/),data_type(ntype))
    irec = irec + 1
end do
end do
end do
print(dzdt(0:(nvar-1),0,9,28,119))
print(dzdt(0,0,9 ,28,99:139:4))
print(dzdt(0,0,10,28,99:139:4))

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
system("rm -f " + fileout(nf))
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "months from " + month_s + " to " + month_e + "daily"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"var","time","lev","lat","lon"/);
dimSizes = (/nvar ,ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"nvar","integer",(/"var"/))
;filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"lev" ,typeof(vars&lev),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"dzdt",typeof(dzdt),(/"var","time","lev","lat","lon"/))
filevardef(ncdf,"dudt",typeof(dudt),(/"var","time","lev","lat","lon"/))
filevardef(ncdf,"dvdt",typeof(dvdt),(/"var","time","lev","lat","lon"/))
filevardef(ncdf,"dTdt",typeof(dTdt),(/"var","time","lev","lat","lon"/))
filevardef(ncdf,"dbdt",typeof(dbdt),(/"var","time","lev","lat","lon"/))


;define var attribute
;filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"lev",vars&lev)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"dzdt",dzdt)
filevarattdef(ncdf,"dudt",dudt)
filevarattdef(ncdf,"dvdt",dvdt)
filevarattdef(ncdf,"dTdt",dTdt)
filevarattdef(ncdf,"dbdt",dbdt)

;ncdf->time = (/vars&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->dzdt = (/dzdt/)
ncdf->dudt = (/dudt/)
ncdf->dvdt = (/dvdt/)
ncdf->dTdt = (/dTdt/)
ncdf->dbdt = (/dbdt/)
end do
end


