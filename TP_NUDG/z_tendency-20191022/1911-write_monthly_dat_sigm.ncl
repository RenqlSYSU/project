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
case_option  = 0 ;0 means NUDG experiment; 1 means noTP experiment
data_type = (/"double","float"/)
ntype = 0
csa = False 

if(case_option.eq.0) then
    case   = (/"CTRL","NUDG6h","NUDG24h"/)
else
    case   = (/"TP_CTRL","TP_CR"/)
end if

;filein  = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_month_Qd_Qeddy_A-g.nc"
;fileout = path + "/project/TP_NUDG/z_tendency-20191022/mdata/" + case + "-Clim_4f6c_month-g600.dat"
filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_month_Qd_Qeddy_A.nc"
fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_Clim_4f6c_month.dat"
;filein  = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_Qd_Qeddy_A.nc"
;fileout = path + "/project/TP_NUDG/observation-20190127/mdata/" + (/"ERA-Interim","NCEP1"/) + "_DJF_4f6c.dat" 
var_name  = (/"f_Qd","f_Qeddy","A"/)  ;this is the variable used to draw, can used as center string
var_name2 = (/"Qd","Qeddy"/)  ;this is the variable used to draw, can used as center string
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4

;area: east Asia
lats = -2  ;-90 ;
latn = 90  ; 90 ;
lonl = 0  ; 0  ;
lonr = 360 ;360 ;
lev2 = (/1000,100/) ;hPa
;lev  = (/600,550,500,450,400,350,300,250,200,150,100/) ;hPa
;lev  = (/1000,925,850,700,600,500,400,300,250,200,150,100/) ;hPa
;lev  = (/1000,975,950,925,900,875,850,825,800,775,750,700,650,600,550,500,450,400,350,300,250,225,200,175,150,125,100/) ;ERA-Interim
lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)
;lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/)
month_s = 1
month_e = 8
nc = 0

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(nc),"r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
;time = cd_calendar(f->time,0)
;select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
;ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)
;nmonth = dimsizes(vars&month)
;forc  = new((/nvar,ntime,nlev,nlat,nlon/),double)
;q3    = new((/3   ,ntime,2   ,nlat,nlon/),double)
ntime   = 4

if(ntype.eq.0) then 
forc_m  = new((/nvar,ntime,nlev,nlat,nlon/),double)
q3_m    = new((/2   ,ntime,2   ,nlat,nlon/),double)
else
forc_m  = new((/nvar,ntime,nlev,nlat,nlon/),float)
q3_m    = new((/2   ,ntime,2   ,nlat,nlon/),float)
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

do nf = nc, nc, 1 
;==========================================================
;read the data
;============================================================
print("handle with " + case(nf))
f      = addfile(filein(nf),"r")
;select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
do nv = 0, nvar-1, 1
;    forc(nv,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    ;forc_m(nv,:,:,:,:) = tofloat(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
    forc_m(nv,:,:,:,:) = f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr})
    ;forc_m(nv,0,:,:,:) = dim_avg_n_Wrap(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}),0)
end do
do nv = 0, nvar-2, 1
;    q3(nv,:,:,:,:)  = f->$var_name2(nv)$(select_time,{(/1000,100/)},{lats:latn},{lonl:lonr})
    ;q3_m(nv,:,:,:,:)  = tofloat(f->$var_name2(nv)$(:,{lev2},{lats:latn},{lonl:lonr}))
    q3_m(nv,:,:,:,:)  = f->$var_name2(nv)$(:,{lev2},{lats:latn},{lonl:lonr})
    ;q3_m(nv,0,:,:,:)  = dim_avg_n_Wrap(f->$var_name2(nv)$(:,{lev2},{lats:latn},{lonl:lonr}),0)
end do
sigm1 := f->sigm1(:,{lev},{lats:latn},{lonl:lonr})
sigm1 := conform(forc_m,sigm1,(/itime,ilev,ilat,ilon/))
;sigm1 := dim_avg_n_Wrap( f->sigm1(:,{lev},{lats:latn},{lonl:lonr}),0)
;forc_m  = where(ismissing(forc_m),0,forc_m)
;q3_m    = where(ismissing(q3_m)  ,0,  q3_m)
print(forc_m(0,0,0,4,4:9))
print(forc_m(1,0,1,4:9,4))

if(csa) then 
forc_m!0 = "nvar"
forc_m!1 = "month"
lev0 = lev
lev := fspan(1000,100,19)
print(lev)
wts = tofloat(lev0) ;weights
wts = 1.0
wts(0) = 0.1 ;1000
wts(1) = 0.2 ;925
wts(2) = 0.5 ;850
knots = nlev - 4
smth = 0.0
nderiv = 0
forc_m_o = csa1x(lev0,forc_m(nvar|:,month|:,lat|:,lon|:,lev|:),wts,knots,smth,nderiv,lev)
 sigm1_o = csa1x(lev0, sigm1(               lat|:,lon|:,lev|:),wts,knots,smth,nderiv,lev)
printVarSummary(forc_m_o)
forc_m_o!0 = "nvar"
forc_m_o!1 = "month"
forc_m_o!2 = "lat"
forc_m_o!3 = "lon"
forc_m_o!4 = "lev"
sigm1_o!0 = "lat"
sigm1_o!1 = "lon"
sigm1_o!2 = "lev"
forc_m := forc_m_o(nvar|:,month|:,lev|:,lat|:,lon|:)
 sigm1 :=  sigm1_o(lev|:,lat|:,lon|:)
nlev = dimsizes(lev)
forc_m@_FillValue = default_fillvalue("float")
 sigm1@_FillValue = default_fillvalue("float")
end if

;linlog = 0 ;1 mean linear interpolation; 0 mean log interpolation
;p_in   = conform(forc_m0,lev0,ilev)
;p_out  = conform(forc_m ,lev ,ilev)
;forc_m = int2p_n(p_in,forc_m0,p_out,linlog,ilev)
;delete([/p_in,p_out/])
;p_in   = conform(sigm10,lev0,ilev-2)
;p_out  = conform(forc_m(0,0,:,:,:) ,lev ,ilev-2)
;sigm1  = int2p_n(p_in,sigm10 ,p_out,linlog,ilev-2)
;delete([/p_in,p_out/])

;ave := wgt_areaave_Wrap(sigm1,1.0,1.0,0)
;ave := dim_avg_n_Wrap(ave,0)
;ave = where(ave.eq.0.0, 1e-5, ave)
;sigm1 := conform(vars(0,:,:,:),ave,ilev-2)
;sigm1 = linmsg_n(sigm1,-1,ilon-1)
;sigm1 = where(sigm1.eq.0.0, 1e-5, sigm1)
;sigm1 := conform(dzdt,sigm1,(/itime,ilev,ilat,ilon/))

;calc ealy(56) and late(78) stage of the summer 
;select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.(month_s+1))
;forc_m(:,0,:,:,:) = dim_avg_n(forc(:,select_time,:,:,:),1)
;q3_m(:,0,:,:,:) = dim_avg_n(q3(:,select_time,:,:,:),1)
;select_time := ind(time(:,1).ge.(month_e-1).and.time(:,1).le.month_e)
;forc_m(:,1,:,:,:) = dim_avg_n(forc(:,select_time,:,:,:),1)
;q3_m(:,1,:,:,:) = dim_avg_n(q3(:,select_time,:,:,:),1)
printVarSummary(forc_m)

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
f0    = todouble(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(lat),ilat))
else
coe110 = tofloat(1/(coslat*coslat*dlon*dlon)) ;111 means z,y,x
coe101 = tofloat(1/dlat/dlat + tanlat/2/dlat)
coe121 = tofloat(1/dlat/dlat - tanlat/2/dlat)
f0    = tofloat(conform(forc_m,2*(2*pi/24.0/3600.0)*sin(lat),ilat))
end if

term1 = center_finite_diff_n(1/sigm1,lev*100,cyclic,opt,ilev)/2/dlev4 
term2 = 1/sigm1/dlev4/dlev4 
coe211 = f0*f0*a*a*(term2 + term1)
coe011 = f0*f0*a*a*(term2 - term1)
coe111 = -2/coslat/coslat/dlon/dlon - 2/dlat/dlat - 2*f0*f0*a*a*term2
coe111 = where(coe111.eq.0,1e-5,coe111)
delete([/term1,term2/])

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
print(coe110(0,0,4,4,4:9))
printVarSummary(coe011)
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

do nv = 0, nvar-2, 1
do nt = 0, ntime-1, 1
do nl = 0, 1, 1
    fbindirwrite(fileout(nf),q3_m(nv,nt,nl,:,:))
end do
end do
end do

do nv = 0, nvar-1, 1
do nt = 0, ntime-1, 1
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
