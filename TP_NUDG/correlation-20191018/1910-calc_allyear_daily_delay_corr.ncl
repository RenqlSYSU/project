;***************************************************************************
;use the post-processing data to calculate JJA daily correlation 
;1. read data PRECL, PRECC, U, V, Q, PS
;2. calc preci, intEKE, intQTE, divg
;3. calc diff between ctrl and nudg
;4. calc JJA daily delay correlation
;5. calc area average and climate average
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/10/18
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
mxlag = 10
fnc = 0

if(case_option.eq.0) then
    case    = (/"CTRL","NUDG6h"/);"NG24h",
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/correlation-20191018/mdata/" + case(fnc) + "-allyear_daily_" + mxlag + "delay_corri.nc"
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_EKE-" + cb + "-" + ca + "day.nc"
end if

var_name= (/"PRECL","PRECC","PS","U","V","Q"/)
g  = 9.8 ;m/(s*s)

area = (/"SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /)
latn = (/ 30  , 45  /)
lonl = (/107  ,115  /)
lonr = (/130  ,130  /)
lev  = (/1000,925,850,700,600,500,400,300,200/)
narea = dimsizes(area)
corr = new((/narea,6,2*mxlag+1/),float)

do na = 0,narea-1,1
print(area(na))
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(3)+".nc","r")
vars:= f->$var_name(3)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = 2
preci := new((/ncase,nyear,ntime,nlat,nlon/),float)
ps    := new((/ncase,nyear,ntime,nlat,nlon/),float)
u     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
q     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
dp    := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
corr1  := new((/6,nlat,nlon,(mxlag+1)/),float)
corr2  := new((/6,nlat,nlon,(mxlag+1)/),float)
corr11 := new((/6,(mxlag+1)/),float)
corr22 := new((/6,(mxlag+1)/),float)

if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nf = 0,1,1
do nt = 0,nyear-1,1
print("read data " + case(nf) + " " + years(nt))
nv = 0
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
preci(nf,nt,:,:,:)   = f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
nv = 1
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
preci(nf,nt,:,:,:)   = preci(nf,nt,:,:,:) + f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
nv = 2
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
ps(nf,nt,:,:,:)   = f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
dp(nf,nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nf,nt,:,:,:),20000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
nv = 3
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
nv = 4
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
nv = 5
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
end do
end do
preci = preci*24*3600*1000 ;convert m/s to mm/day
printVarSummary(u)

div850 := uv2dv_cfd( u(:,:,:,2,:,:), v(:,:,:,2,:,:), vars&lat, vars&lon, 2)
if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb    
    opt = False
    edd_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    edd_v = bw_bandpass_filter(v,fca,fcb,opt,1)
    edd_q = bw_bandpass_filter(q,fca,fcb,opt,1)
else
    edd_q = q - conform(q,dim_avg_n(q,1),(/0,2,3,4,5/))
    edd_u = u - conform(u,dim_avg_n(u,1),(/0,2,3,4,5/))
    edd_v = v - conform(v,dim_avg_n(v,1),(/0,2,3,4,5/))   
    preci = preci  - conform(preci ,dim_avg_n(preci ,1),(/0,2,3,4/))   
    div850= div850 - conform(div850,dim_avg_n(div850,1),(/0,2,3,4/))   
end if
intEKE := dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,3)/g/2.0   ;ps-200hPa integral get (nyear,ntime,nlat,nlon)
QTE := uv2dv_cfd( edd_q*edd_u, edd_q*edd_v, vars&lat, vars&lon, 2)
intQTE := dim_sum_n(QTE*dp,3)/g 
delete([/edd_u,edd_v,edd_q/])

preci(1,:,:,:,:)  = preci(1,:,:,:,:)  - preci(0,:,:,:,:)
intEKE(1,:,:,:,:) = intEKE(1,:,:,:,:) - intEKE(0,:,:,:,:)
intQTE(1,:,:,:,:) = intQTE(1,:,:,:,:) - intQTE(0,:,:,:,:)
div850(1,:,:,:,:) = div850(1,:,:,:,:) - div850(0,:,:,:,:)

preci  := reshape(preci, (/ncase,nyear*ntime,nlat,nlon/))
intEKE := reshape(intEKE,(/ncase,nyear*ntime,nlat,nlon/))
intQTE := reshape(intQTE,(/ncase,nyear*ntime,nlat,nlon/))
div850 := reshape(div850,(/ncase,nyear*ntime,nlat,nlon/))

copy_VarMeta(vars(0,0,:,:),preci(0,0,:,:))
preci!0 = "case"
preci!1 = "time"
copy_VarMeta(preci,intEKE)
copy_VarMeta(preci,intQTE)
copy_VarMeta(preci,div850)
corr1(0,:,:,:) = esccr(preci(case|fnc,lat|:,lon|:,time|:),intEKE(case|fnc,lat|:,lon|:,time|:),mxlag)
corr1(1,:,:,:) = esccr(preci(case|fnc,lat|:,lon|:,time|:),intQTE(case|fnc,lat|:,lon|:,time|:),mxlag)
corr1(2,:,:,:) = esccr(preci(case|fnc,lat|:,lon|:,time|:),div850(case|fnc,lat|:,lon|:,time|:),mxlag)
corr1(3,:,:,:) = esccr(intQTE(case|fnc,lat|:,lon|:,time|:),intEKE(case|fnc,lat|:,lon|:,time|:),mxlag)
corr1(4,:,:,:) = esccr(div850(case|fnc,lat|:,lon|:,time|:),intEKE(case|fnc,lat|:,lon|:,time|:),mxlag)
corr1(5,:,:,:) = esccr(div850(case|fnc,lat|:,lon|:,time|:),intQTE(case|fnc,lat|:,lon|:,time|:),mxlag)

corr2(0,:,:,:) = esccr(intEKE(case|fnc,lat|:,lon|:,time|:),preci(case|fnc,lat|:,lon|:,time|:),mxlag)
corr2(1,:,:,:) = esccr(intQTE(case|fnc,lat|:,lon|:,time|:),preci(case|fnc,lat|:,lon|:,time|:),mxlag)
corr2(2,:,:,:) = esccr(div850(case|fnc,lat|:,lon|:,time|:),preci(case|fnc,lat|:,lon|:,time|:),mxlag)
corr2(3,:,:,:) = esccr(intEKE(case|fnc,lat|:,lon|:,time|:),intQTE(case|fnc,lat|:,lon|:,time|:),mxlag)
corr2(4,:,:,:) = esccr(intEKE(case|fnc,lat|:,lon|:,time|:),div850(case|fnc,lat|:,lon|:,time|:),mxlag)
corr2(5,:,:,:) = esccr(intQTE(case|fnc,lat|:,lon|:,time|:),div850(case|fnc,lat|:,lon|:,time|:),mxlag)

copy_VarMeta(vars(0,0,:,:),corr1(0,:,:,0))
copy_VarMeta(vars(0,0,:,:),corr2(0,:,:,0))
do nm = 0,mxlag,1
    corr11(:,nm) = wgt_areaave(corr1(:,:,:,nm),1.0,1.0,0) ;The rightmost dimensions should correspond to "latitude" (lat) and "longitude" (lon)
    corr22(:,nm) = wgt_areaave(corr2(:,:,:,nm),1.0,1.0,0)
end do
corr(na,:,0:mxlag-1) = corr22(:,1:mxlag:-1) ;negative lag
corr(na,:,mxlag:)    = corr11(:,0:mxlag:)   ;positive lag
end do

;=======================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "6 delay correlation between intEKE and preci, intQTE and preci, div850 and preci" + \
                        "intQTE and intEKE, div850 and intEKE, div850 and intQTE"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"area","nmber","delay"/);
dimSizes = (/narea,6    ,2*mxlag+1/) ;
dimUnlim = (/False,False,False      /);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"area" ,"integer",(/"area" /))
filevardef(ncdf,"nmber","integer",(/"nmber"/))
filevardef(ncdf,"delay","integer",(/"delay"/))
filevardef(ncdf,"corr" ,"float"  ,(/"area","nmber","delay"/)) 

ncdf->narea  = ispan(1,2,1)
ncdf->nmber  = ispan(1,6,1)
ncdf->delay  = ispan(-mxlag,mxlag,1)
ncdf->corr   = (/corr/)
end

