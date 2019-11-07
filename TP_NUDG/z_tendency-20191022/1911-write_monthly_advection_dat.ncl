;***************************************************************************
;use the 3 advection forcing to calc dzdt by SOR (Successive Over Relaxation)
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

filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_MJJA_daily_adevection.nc"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_MJJA_3advection_monthly.dat"
var_name = (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4

;area: east Asia
lats = -2 ;-90 ;
latn = 90  ; 90 ;
lonl = 0  ; 0  ;
lonr = 360 ;360 ;
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
month_s = 5
month_e = 8

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0),"r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)
forc  = new((/nvar,ntime,nlev,nlat,nlon/),double)
forc_m= new((/nvar,2    ,nlev,nlat,nlon/),double)

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
dlev3  = conform(vars,dlev,ilev-2)
coslat = conform(vars,cos(lat),ilat-2)
tanlat = conform(vars,tan(lat),ilat-2)
coe110 = 1/(coslat*coslat*dlon*dlon) ;111 means z,y,x
coe101 = 1/dlat/dlat + tanlat/2/dlat
coe121 = 1/dlat/dlat - tanlat/2/dlat

opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 1, 1, 1 
;==========================================================
;read the data
;============================================================
print("handle with " + case(nf))
f      = addfile(filein(nf),"r")
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
do nv = 0, nvar-1, 1
    forc(nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
sigm0 = f->sigm0(select_time,{lev},{lats:latn},{lonl:lonr})
forc  = linmsg_n(forc,-1,ilon) ;-1 means the beginning and end points will be set to the nearest non-missing value.

ave := wgt_areaave_Wrap(sigm0,1.0,1.0,0)
ave := dim_avg_n_Wrap(ave,0)
ave  = where(ave.eq.0.0, 1e-12, ave)
sigm1 = conform(vars,ave*R/lev/100,ilev-2)

;calc ealy(56) and late(78) stage of the summer 
select_time := ind(time(:,1).ge.5.and.time(:,1).le.6)
forc_m(:,0,:,:,:) = dim_avg_n(forc(:,select_time,:,:,:),1)
select_time := ind(time(:,1).ge.7.and.time(:,1).le.8)
forc_m(:,1,:,:,:) = dim_avg_n(forc(:,select_time,:,:,:),1)
printVarSummary(forc_m)

if (any(ismissing(forc))) then
    print("forc has missing")
end if
if (any(ismissing(sigm1))) then
    print("sigm1 has missing")
end if

;==========================================================
;define coefficient and calc dzdt
;============================================================
f0    = conform(vars,2*(2*pi/24.0/3600.0)*sin(lat),ilat-2)
term1 = center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev-2)/2/dlev3 
term2 = 1/sigm1/dlev3/dlev3 
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

print(forc(0,1,4,4,4:9))
print(coe110(4,4,4:9))

;==========================================================
;write to binary file
;============================================================
setfileoption("bin","WriteByteOrder","LittleEndian")
system("rm -f " + fileout(nf))

do nv = 0, nvar-1, 1
do nt = 0, 1, 1
do nl = 0, nlev-1, 1
    fbindirwrite(fileout(nf),forc_m(nv,nt,nl,:,:))
end do
end do
end do

do nl = 0, nlev-1, 1
    fbindirwrite(fileout(nf),coe111(nl,:,:))
    fbindirwrite(fileout(nf),coe110(nl,:,:))
    fbindirwrite(fileout(nf),coe101(nl,:,:))
    fbindirwrite(fileout(nf),coe121(nl,:,:))
    fbindirwrite(fileout(nf),coe011(nl,:,:))
    fbindirwrite(fileout(nf),coe211(nl,:,:))
    fbindirwrite(fileout(nf),f0(nl,:,:))
end do

end do
end
