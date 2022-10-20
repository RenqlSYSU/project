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
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case   = (/"CTRL","NUDG6h","NUDG24h"/)
else
    case   = (/"TP_CTRL","TP_CR"/)
end if

filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_Qd_Qeddy_A.nc"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_dzdt.nc"
var_name  = (/"f_Qd","f_Qd_t","f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
var_name2 = (/"Qd","Qd_t","Qeddy"/)  ;this is the variable used to draw, can used as center string
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats = -20 ;-90 ;
latn = 60  ; 90 ;
lonl = 20  ; 0  ;
lonr = 220 ;360 ;
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
month_s = 5
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0),"r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)
forc  = new((/nvar,ntime,nlev,nlat,nlon/),double)
q3    = new((/3   ,ntime,2   ,nlat,nlon/),double)
dzdt  = new((/nvar,ntime,nlev,nlat,nlon/),double)
term  = new((/nvar,ntime,nlev,nlat,nlon/),double)
diff  = new((/nvar/),double)

rf   = 0.25;0.75 ;relaxing factor
critical = 1e-8
itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0

dlon = lon(4)-lon(3)
dlat = lat(4)-lat(3)
dlev = lev
dlev(1:(nlev-1)) = (- lev(0:(nlev-2)) + lev(1:(nlev-1)))*100 
dlev(0) = dlev(1) 
dlev5  = conform(dzdt,dlev,ilev)
coslat = conform(dzdt,cos(lat),ilat)
tanlat = conform(dzdt,tan(lat),ilat)
coe110 = 1/(coslat*coslat*dlon*dlon) ;111 means z,y,x
coe101 = 1/dlat/dlat + tanlat/2/dlat
coe121 = 1/dlat/dlat - tanlat/2/dlat

opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 0, ncase-1, 1 
;==========================================================
;read the data
;============================================================
print("handle with " + case(nf))
f      = addfile(filein(nf),"r")
do nv = 0, nvar-1, 1
    forc(nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
do nv = 0, 2, 1
    q3(nv,:,:,:,:)  = f->$var_name2(nv)$(select_time,{(/1000,100/)},{lats:latn},{lonl:lonr})
end do
sigm1 = f->sigm1(select_time,{lev},{lats:latn},{lonl:lonr})
copy_VarMeta(forc,dzdt)
;forc  = where(ismissing(forc),0,forc)
forc = linmsg_n(forc,-1,ilon) ;-1 means the beginning and end points will be set to the nearest non-missing value.

ave = wgt_areaave_Wrap(sigm1,1.0,1.0,0)
ave := dim_avg_n_Wrap(ave,0)
ave = where(ave.eq.0.0, 1e-5, ave)
sigm1 := conform(dzdt,ave,ilev)
if (any(ismissing(ave))) then
    print("ave has missing")
end if
;sigm1 = linmsg_n(sigm1,-1,ilon-1)
;sigm1 = where(sigm1.eq.0.0, 1e-5, sigm1)
;sigm1 := conform(dzdt,sigm1,(/itime,ilev,ilat,ilon/))

if (any(ismissing(forc))) then
    print("forc has missing")
end if
if (any(ismissing(sigm1))) then
    print("sigm1 has missing")
end if
printVarSummary(forc)
printVarSummary(sigm1)

;==========================================================
;define coefficient and calc dzdt
;============================================================
f0    = conform(dzdt,2*(2*pi/24.0/3600.0)*sin(lat),ilat)
term1 = center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev)/2/dlev5 
term2 = 1/sigm1/dlev5/dlev5 
print(term1(0,0,0,0,0:9))
print(term2(0,0,0,0,0:9))
coe211 = f0*f0*a*a*(term2 + term1)
coe011 = f0*f0*a*a*(term2 - term1)
coe111 = -2/coslat/coslat/dlon/dlon - 2/dlat/dlat - 2*f0*f0*a*a*term2
coe111 = where(coe111.eq.0,1e-5,coe111)
delete([/term1,term2/])

if (any(ismissing(coe111))) then
    print("coe111 has missing")
end if
if (any(ismissing(coe110))) then
    print("coe110 has missing")
end if
if (any(ismissing(coe101))) then
    print("coe101 has missing")
end if
if (any(ismissing(coe121))) then
    print("coe121 has missing")
end if
if (any(ismissing(coe011))) then
    print("coe211 has missing")
end if
if (any(ismissing(coe211))) then
    print("coe211 has missing")
end if

dzdt = 0.0
iter = 1
printVarSummary(dzdt)
if (any(ismissing(dzdt))) then
    print("dzdt has missing")
end if

do while (iter.lt.5000)
    print("")
    print(iter + " iteration")
    dzdt0 = dzdt
    do nz = 1, nlev-2, 1
    do ny = 1, nlat-2, 1
    do nx = 1, nlon-2, 1
        term(:,:,nz,ny,nx) = coe111(:,:,nz,ny,nx)*dzdt(:,:,nz,ny,nx) + \
                             coe110(:,:,nz,ny,nx)*dzdt(:,:,nz,ny,nx+1) + coe121(:,:,nz,ny,nx)*dzdt(:,:,nz,ny+1,nx) + coe211(:,:,nz,ny,nx)*dzdt(:,:,nz+1,ny,nx) +\
                             coe110(:,:,nz,ny,nx)*dzdt(:,:,nz,ny,nx-1) + coe101(:,:,nz,ny,nx)*dzdt(:,:,nz,ny-1,nx) + coe011(:,:,nz,ny,nx)*dzdt(:,:,nz-1,ny,nx) 
        dzdt(:,:,nz,ny,nx) = dzdt(:,:,nz,ny,nx) + (rf/coe111(:,:,nz,ny,nx))*(f0(:,:,nz,ny,nx)*a*a*forc(:,:,nz,ny,nx)-term(:,:,nz,ny,nx))
        ;print("")
    end do
    end do
    end do
    nx = 30
    ny = 10
    nz = 10
    print(term(2,2,nz,ny,nx)+" = "+coe110(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz,ny,nx+1)+" + "+coe121(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz,ny+1,nx)+" + "+coe211(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz+1,ny,nx) +\
                             " + "+coe110(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz,ny,nx-1)+" + "+coe101(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz,ny-1,nx)+" + "+coe011(2,2,nz,ny,nx)+"*"+dzdt(2,2,nz-1,ny,nx))
    print(dzdt(2,2,nz,ny,nx)+" = "+dzdt(2,2,nz,ny,nx)+" + ( "+rf+" / "+coe111(2,2,nz,ny,nx)+" ) * ("+f0(2,2,nz,ny,nx)+"*a*a*"+forc(2,2,nz,ny,nx)+" - "+term(2,2,nz,ny,nx)+" )")

;boundary conditions
    dzdt(:,:,:,:,0     ) = dzdt(:,:,:,:,1     )
    dzdt(:,:,:,:,nlon-1) = dzdt(:,:,:,:,nlon-2)
    
    dzdt(:,:,:,0     ,:) = 0 ;dzdt(:,:,:,1     ,:)
    dzdt(:,:,:,nlat-1,:) = dzdt(:,:,:,nlat-2,:)
    
    dzdt(3,:,0     ,:,:) = dzdt(3,:,1     ,:,:)  ;lower boundary for A
    dzdt(3,:,nlev-1,:,:) = dzdt(3,:,nlev-2,:,:)  ;upper boundary for A
    dzdt(0:2,:,0     ,:,:) = dzdt(0:2,:,1     ,:,:) + (R*dlev(0)/lev(0)/100)*forc(0:2,:,0,:,:)    ;lower boundary for Qd, Qd_t,Qeddy
    dzdt(0:2,:,nlev-1,:,:) = dzdt(0:2,:,nlev-2,:,:) - (R*dlev(nlev-1)/lev(nlev-1)/100)*forc(0:2,:,nlev-1,:,:) ;upper boundary for Qd, Qd_t,Qeddy
    
if (all(ismissing(term))) then
    print("term is all missing")
end if
    
    ;printVarSummary(dzdt)
    ;printVarSummary(term)
    do nv = 0, nvar-1 ,1 
    if (all(ismissing(dzdt(nv,:,:,:,:)))) then
        print("dzdt induced by " + var_name(nv) + " is all missing")
    end if
    diff1 = abs(dzdt(nv,:,:,:,:)-dzdt0(nv,:,:,:,:))
    if (all(ismissing(diff1))) then
        print("diff of " + var_name(nv) + " is all missing")
    end if
    copy_VarMeta(forc(0,:,:,:,:),diff1)
    diff(nv) = max(diff1)
    print("Iterative difference of " + var_name(nv) + " " + diff(nv))
    end do
    if (max(diff).lt.critical) then 
        break
    end if 
    iter = iter + 1
end do

do ilat = 0, nlat-1
if (abs(lat(ilat)*180.0/pi ).lt. 10. ) then
f0(:,:,ilat,:)= f0@_FillValue
end if
end do

dudt = - center_finite_diff_n(dzdt,lat    ,cyclic,opt,ilat)/f0
dvdt =   center_finite_diff_n(dzdt,lon    ,cyclic,opt,ilon)/f0
dTdt = - center_finite_diff_n(dzdt,lev*100,cyclic,opt,ilev)*conform(dzdt,lev*100,ilev)/R
dbdt = - center_finite_diff_n(dTdt,lat    ,cyclic,opt,ilat)

copy_VarMeta(forc,dzdt)
copy_VarMeta(forc,dudt)
copy_VarMeta(forc,dvdt)
copy_VarMeta(forc,dTdt)
copy_VarMeta(forc,dbdt)

nvar@long_name = "dzdt induced by f_Qd, f_Qd_t, f_Qeddy and A" 
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
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "months from " + month_s + " to " + month_e + "daily"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"nvar","time","lev","lat","lon"/);
dimSizes = (/nvar  ,ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(u&time),getvardims(u&time))
filevardef(ncdf,"lev" ,typeof(u&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
filevardef(ncdf,"dzdt",typeof(dzdt),(/"nvar","time","lev","lat","lon"/))
filevardef(ncdf,"dudt",typeof(dudt),(/"nvar","time","lev","lat","lon"/))
filevardef(ncdf,"dvdt",typeof(dvdt),(/"nvar","time","lev","lat","lon"/))
filevardef(ncdf,"dTdt",typeof(dTdt),(/"nvar","time","lev","lat","lon"/))
filevardef(ncdf,"dbdt",typeof(dbdt),(/"nvar","time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"nvar",nvar)
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lev",u&level)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"dzdt",dzdt)
filevarattdef(ncdf,"dudt",dudt)
filevarattdef(ncdf,"dvdt",dvdt)
filevarattdef(ncdf,"dTdt",dTdt)
filevarattdef(ncdf,"dbdt",dbdt)

ncdf->nvar = ispan(1,nvar,1) 
ncdf->time = (/u&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->dzdt = (/dzdt/)
ncdf->dudt = (/dudt/)
ncdf->dvdt = (/dvdt/)
ncdf->dTdt = (/dTdt/)
ncdf->dbdt = (/dbdt/)
end do
end


