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

year_option = True ;False ;True means yearly
ave_sigm = False ;recalculate the f_Qd and f_Qeddy use fixed f
nyear = 37
decomp = False;True ;

case   = "ERA-Interim"
;case   = "NCEP1"

if(.not.(decomp)) then 
filein  = path + case + "_DJF_forc_re-year.nc"
fileout = path + case + "_DJF_3b6c.dat"
var_name  = (/"turn","b1z","b4"/)  ;this is the variable used to draw, can used as center string
;var_name  = (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
else
filein  = path + case + "_DJF_forc_in2-year.nc"
fileout = path + case + "_DJF_3b6c.dat"
var_name  = (/"b3x","b3y","b1x","b1y"/)  ;this is the variable used to draw, can used as center string
var_name2 = (/"t3x","t3y"/)  
end if

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4

;area: east Asia
lats = 10  ;-90 ;
latn = 90  ; 90 ;
lonl = 0   ; 30  ; 0  ;
lonr = 360 ;150 ;360 ;
lev  = (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,175,150,125, 100/)
lev2 = (/1000,100/) ;hPa

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
forc_m  = new((/nvar,nyear,nlev,nlat,nlon/),double)
else
forc_m  = new((/nvar,nyear,nlev,nlat,nlon/),float)
end if

itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0

opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data
;============================================================
f  = addfile(filein,"r")
do nv = 0, nvar-1, 1
if(ntype.eq.0) then 
    forc_m(nv,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
else
    forc_m(nv,:,:,:,:) = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
end if
end do

if(nv.eq.4) then 
q3_m    = new((/2   ,nyear,2   ,nlat,nlon/),double)
do nv = 0, 1, 1
    q3_m(nv,:,:,:,:)  = f->$var_name2(nv)$(:,{lev2},{lats:latn},{lonl:lonr})
end do
q3_m = -q3_m
end if

f11 = addfile(path + case + "_DJF_forc_in-year.nc","r")
if(ntype.eq.0) then 
    sigm0 := dim_avg_n_Wrap(f11->sigm0(:,{lev},{lats:latn},{lonl:lonr}),0)
else
    sigm0 := tofloat(dim_avg_n_Wrap(f11->sigm0(:,{lev},{lats:latn},{lonl:lonr}),0))
end if
sigm0 := conform(forc_m,wgt_areaave(sigm0, 1.0, 1.0, 0),(/ilev/))
sigm1 = sigm0*R/conform(forc_m,lev*100,ilev)
print(forc_m(0,0,0,4,4:9))
print(forc_m(1,0,1,4:9,4))

if (any(ismissing(forc_m))) then
    print("forc has missing")
end if
if (any(ismissing(sigm1))) then
    print("sigm1 has missing")
end if

;==========================================================
;define coefficient and calc dzdt
;============================================================
dlon = lon(4)-lon(3)
dlat = lat(4)-lat(3)
dlev = lev
dlev(1:(nlev-1)) = (- lev(0:(nlev-2)) + lev(1:(nlev-1)))*100 
dlev(0) = dlev(1) 
dlev4  = conform(forc_m,dlev,ilev)
coslat = conform(forc_m,cos(lat),ilat)
tanlat = conform(forc_m,tan(lat),ilat)
if(ntype.eq.0) then 
    coe110 = todouble(1/(coslat*coslat*dlon*dlon)) ;111 means z,y,x
    coe101 = todouble(1/dlat/dlat + tanlat/2/dlat)
    coe121 = todouble(1/dlat/dlat - tanlat/2/dlat)
    if(ave_sigm) then
    f0    = todouble(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(45.0*pi/180.0),ilat));45.0*pi/180.0
    else
    f0    = todouble(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(lat),ilat));45.0*pi/180.0
    end if
    term1 = center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev)/2/dlev4 
    term2 = 1/sigm1/dlev4/dlev4 
    coe211 = f0*f0*a*a*(term2 + term1);
    coe011 = f0*f0*a*a*(term2 - term1);
    coe111 = -2/coslat/coslat/dlon/dlon - 2/dlat/dlat - 2*f0*f0*a*a*term2
else
    coe110 = tofloat(1/(coslat*coslat*dlon*dlon)) ;111 means z,y,x
    coe101 = tofloat(1/dlat/dlat + tanlat/2/dlat)
    coe121 = tofloat(1/dlat/dlat - tanlat/2/dlat)
    f0    = tofloat(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(lat),ilat))
    term1 = tofloat(center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev)/2/dlev4 )
    term2 = tofloat(1/sigm1/dlev4/dlev4) 
    coe211 = tofloat(f0*f0*a*a*(term2 + term1));
    coe011 = tofloat(f0*f0*a*a*(term2 - term1));
    coe111 = tofloat(-2/coslat/coslat/dlon/dlon - 2/dlat/dlat - 2*f0*f0*a*a*term2)
end if

coe111 = where(coe111.eq.0,1e-5,coe111)
delete([/term1,term2,dlev4/])

coe111@_FillValue = default_fillvalue(data_type(ntype))
coe110@_FillValue = default_fillvalue(data_type(ntype))
coe101@_FillValue = default_fillvalue(data_type(ntype))
coe121@_FillValue = default_fillvalue(data_type(ntype))
coe011@_FillValue = default_fillvalue(data_type(ntype))
coe211@_FillValue = default_fillvalue(data_type(ntype))
f0@_FillValue = default_fillvalue(data_type(ntype))

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

print(forc_m(0,0,0,4,4:9))
print(forc_m(1,0,1,4:9,4))
print(coe111(0,0,4,4,4:9))
printVarSummary(coe211)
printVarSummary(f0)

;==========================================================
;write to binary file
;============================================================
setfileoption("bin","WriteByteOrder","LittleEndian")
system("rm -f " + fileout)

do nv = 0, nvar-1, 1
do nt = 0, nyear-1, 1
do nl = 0, nlev-1, 1
    fbindirwrite(fileout,forc_m(nv,nt,nl,:,:))
end do
end do
end do

if(nv.eq.4) then 
do nv = 0, 1, 1
do nt = 0, nyear-1, 1
do nl = 0, 1, 1
    fbindirwrite(fileout,q3_m(nv,nt,nl,:,:))
end do
end do
end do
end if

do nv = 0, nvar-1, 1  ;
do nt = 0, nyear-1, 1 ;
do nl = 0, nlev-1, 1
    fbindirwrite(fileout,coe111(nv,nt,nl,:,:))
    fbindirwrite(fileout,coe110(nv,nt,nl,:,:))
    fbindirwrite(fileout,coe101(nv,nt,nl,:,:))
    fbindirwrite(fileout,coe121(nv,nt,nl,:,:))
    fbindirwrite(fileout,coe011(nv,nt,nl,:,:))
    fbindirwrite(fileout,coe211(nv,nt,nl,:,:))
    fbindirwrite(fileout,f0(nv,nt,nl,:,:))
end do
end do
end do
end
