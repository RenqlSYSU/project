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
data_type = (/"double","float"/)
ntype = 1

obs_option = False ;True
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
year_option = False;True
csa = False 
ave_sigm = False;True 
nc = 1
ntime = 4
month_s = 1
month_e = 8

if(case_option.eq.0) then
    case   = (/"CTRL","NUDG6h","NUDG24h"/)
end if

if(case_option.eq.1) then
    case   = (/"TP_CTRL","TP_CR"/)
end if

if(case_option.eq.2) then
    case    = (/"F2000","NUDG","NG58"/)
end if

if(.not.(obs_option)) then 
if(ntime.gt.12) then 
filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_daily_Qd_Qeddy_A.nc"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_4f6c_month.dat"
end if
if(ntime.eq.4) then 
filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_advect2.nc"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_3b6c_month.dat"
end if
if(ntime.eq.1) then 
filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_DJF_Qd_Qeddy_A.nc"
fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_4f6c_month.dat"
end if
else
if(ntime.eq.4) then 
filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_month_Qd_Qeddy_A-geo.nc"
fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_4f6c_month.dat"
else
filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_Qd_Qeddy_A_1204.nc"
fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_4f6c.dat" 
end if
end if
var_name  = (/"b1","b2","b3"/)  ;this is the variable used to draw, can used as center string
ncase   = dimsizes(case)
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
;lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)
lev  = (/1000,950,900,850,800, 750,700,650,600,550, 500,450,400,350,300, 250,200,150,100/) ;19levels, hpa

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
print(filein(nc))
f    = addfile(filein(nc),"r")
if(ntime.eq.1) then 
vars:= f->$var_name(0)$({lev},{lats:latn},{lonl:lonr})
else
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
end if
;time = cd_calendar(f->time,0)
;select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
;ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)

if(ntype.eq.0) then 
forc_m  = new((/nvar,ntime,nlev,nlat,nlon/),double)
else
forc_m  = new((/nvar,ntime,nlev,nlat,nlon/),float)
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

do nf = 0, nc, 1
;==========================================================
;read the data
;============================================================
print("handle with " + case(nf))
f      = addfile(filein(nf),"r")
;select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
do nv = 0, nvar-1, 1
if(year_option) then 
    forc_m(nv,0,:,:,:) = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}),0)
else
;    forc(nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
if(ntime.eq.1) then 
    forc_m(nv,0,:,:,:) = f->$var_name(nv)$({lev},{lats:latn},{lonl:lonr})
else
    if(ntype.eq.0) then 
    forc_m(nv,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    else
    forc_m(nv,:,:,:,:) = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    end if
end if
end if
end do

if(year_option) then 
sigm0 := dim_avg_n_Wrap( f->sigm0(:,{lev},{lats:latn},{lonl:lonr}),0)
sigm0 := conform(forc_m,wgt_areaave(sigm0, 1.0, 1.0, 0),(/ilev/))
else
    if(ntime.eq.1) then 
    sigm0 := f->sigm0({lev},{lats:latn},{lonl:lonr})
    sigm0 := conform(forc_m,wgt_areaave(sigm0, 1.0, 1.0, 0),(/ilev/))
    end if
    if(ntime.eq.4) then 
    sigm0 := f->sigm0(:,{lev},{lats:latn},{lonl:lonr})
    sigm0 := conform(forc_m,wgt_areaave(sigm0, 1.0, 1.0, 0),(/itime,ilev/))
    end if
    if(ntime.gt.12) then 
        if(ntype.eq.0) then 
        sigm0 := dim_avg_n_Wrap(f->sigm0(0:4,{lev},{lats:latn},{lonl:lonr}),0)
        else
        sigm0 := tofloat(dim_avg_n_Wrap(f->sigm0(0:4,{lev},{lats:latn},{lonl:lonr}),0))
        end if
    sigm0 := conform(forc_m,wgt_areaave(sigm0, 1.0, 1.0, 0),(/ilev/))
    end if
end if
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
    f0     = todouble(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(lat),ilat));45.0*pi/180.0
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
system("rm -f " + fileout(nf))

do nv = 0, nvar-1, 1
do nt = 0, ntime-1, 1
do nl = 0, nlev-1, 1
    fbindirwrite(fileout(nf),forc_m(nv,nt,nl,:,:))
end do
end do
end do

do nv = 0, nvar-1, 1  ;
do nt = 0, ntime-1, 1 ;
do nl = 0, nlev-1, 1
    fbindirwrite(fileout(nf),coe111(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),coe110(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),coe101(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),coe121(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),coe011(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),coe211(nv,nt,nl,:,:))
    fbindirwrite(fileout(nf),f0(nv,nt,nl,:,:))
end do
end do
end do

end do
end
