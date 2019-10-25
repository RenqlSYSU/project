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
var_name = (/"f_Qd","f_Qd_t","f_Qeddy","A","sigm1"/)  ;this is the variable used to draw, can used as center string
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
dzdt  = new((/nvar,ntime,nlev,nlat,nlon/),double)
term  = new((/nvar,ntime,nlev,nlat,nlon/),double)
diff  = new((/nvar/),double)

itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
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
do nv = 0, nvar-2, 1
    forc(nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
nv = nvar - 1
sigm1 = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
printVarSummary(w)

;==========================================================
;define coefficient and calc dzdt
;============================================================
rf   = 1.25 ;relaxing factor
lat  = lat*pi/180.0
lon  = lon*pi/180.0
dlon = lon(4)-lon(3)
dlat = lat(4)-lat(3)
dlev = lev
dlev(1:(nlev-1)) = (lev(0:(nlev-2)) - lev(1:(nlev-1)))*100 
dlev(0) = dlev(1) 
f0     = conform(vars(0,:,:,:), 2*(2*pi/24.0/3600.0)*sin(lat), ilat-2)
coslat = conform(vars(0,:,:,:),cos(lat),ilat-2)
tanlat = conform(vars(0,:,:,:),tan(lat),ilat-2)
dlev5  = conform(dzdt,dlev,ilev)
sigm1 := conform(dzdt,sigm1,(/itime,ilev,ilat,ilon/))

coe110 = 1/(coslat*coslat*dlon*dlon) ;111 means z,y,x
coe101 = 1/dlat/dlat + tanlat/2/dlat
coe121 = 1/dlat/dlat - tanlat/2/dlat
term1 = center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev)/2/dlev5 
term2 = 1/sigm1/dlev5/dlev5 
coe211 = f0*f0*a*a*(term2 + term1)
coe011 = f0*f0*a*a*(term2 - term1)
coe111 = -2/coslat/coslat/dlon/dlon - 2/dlat/dlat - 2*f0*f0*a*a*term2
delete([/term1,term2/])

dzdt = 1
iter = 1
do while (iter.lt.200)
    print("")
    print(iter + "iteration")
    dzdt0 = dzdt
    do nz = 1, nlev-2, 1
    do ny = 1, nlat-2, 1
    do nx = 1, nlon-2, 1
        term(:,:,nz,ny,nx) = coe110(nz,ny,nx)*dzdt(:,:,nz,ny,nx+1) + coe121(nz,ny,nx)*dzdt(:,:,nz,ny+1,nx) + coe211(nz,ny,nx)*dzdt(:,:,nz+1,ny,nx) +\
                             coe110(nz,ny,nx)*dzdt(:,:,nz,ny,nx-1) + coe101(nz,ny,nx)*dzdt(:,:,nz,ny-1,nx) + coe011(nz,ny,nx)*dzdt(:,:,nz-1,ny,nx) 
        dzdt(:,:,nz,ny,nx) = dzdt(:,:,nz,ny,nx) + (rf/coe111(nz,ny,nx))*(f0(nz,ny,nx)*a*a*forc(:,:,nz,ny,nx)-term(:,:,nz,ny,nx))
    end do
    end do
    end do
    
    do nv = 0, nvar-1 ,1 
        diff(nv) = max(abs(dzdt(nv,:,:,:,:)-dzdt0(nv,:,:,:,:)))
        print("Iterative difference of " + var_name(nv) + " " + diff(nv))
    end do
    if (max(diff).lt.critical) then 
        break
    end if 
    k = k+1   
end do
;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "months from " + month_s + " to " + month_e + "daily three forcing factors, Qd, Qeddy, A"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lev","lat","lon"/);
dimSizes = (/ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(u&time),getvardims(u&time))
filevardef(ncdf,"lev" ,typeof(u&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(u&lat),getvardims(u&lat))
filevardef(ncdf,"lon" ,typeof(u&lon),getvardims(u&lon))
filevardef(ncdf,"Qeddy",typeof(Qeddy),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qd"   ,typeof(Qd),(/"time","lev","lat","lon"/))
filevardef(ncdf,"Qd_t" ,typeof(Qd_t),(/"time","lev","lat","lon"/))
filevardef(ncdf,"A"    ,typeof(A),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qd",typeof(f_Qd),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qd_t",typeof(f_Qd_t),(/"time","lev","lat","lon"/))
filevardef(ncdf,"f_Qeddy",typeof(f_Qeddy),(/"time","lev","lat","lon"/))


;define var attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lev",u&level)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"Qeddy",Qeddy)
filevarattdef(ncdf,"Qd"   ,Qd)
filevarattdef(ncdf,"Qd_t" ,Qd_t)
filevarattdef(ncdf,"A"    ,A)
filevarattdef(ncdf,"f_Qeddy",f_Qeddy)
filevarattdef(ncdf,"f_Qd"   ,f_Qd)
filevarattdef(ncdf,"f_Qd_t" ,f_Qd_t)

ncdf->time = (/u&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->Qeddy = (/Qeddy/)
ncdf->Qd    = (/Qd/)
ncdf->Qd_t  = (/Qd_t/)
ncdf->A     = (/A/)
ncdf->f_Qeddy = (/f_Qeddy/)
ncdf->f_Qd    = (/f_Qd/)
ncdf->f_Qd_t  = (/f_Qd_t/)
end do
end


