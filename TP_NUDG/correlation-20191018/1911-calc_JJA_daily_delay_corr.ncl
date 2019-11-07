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
month_s = 3
month_e = 7

if(case_option.eq.0) then
    case    = (/"CTRL","NUDG6h"/);"NG24h",
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/correlation-20191018/mdata/" + case(fnc) + "-daily" + mxlag + "delay_corri_" + month_s + "-" + month_e + ".nc"
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_month_ave_EKE-" + cb + "-" + ca + "day.nc"
end if

var_name= (/"PRECL","PRECC","U","V","Q","Q1"/)
varname = (/"preci","QTE850","EKE850","EKE200","q1_850","q1_200","div850","div200"/)
nvar2 = dimsizes(varname)
nvar3 = (nvar2-1)*nvar2/2
print(nvar2)
print(nvar3)
g  = 9.8 ;m/(s*s)

area = (/"SEA","NEA"/); "TP","STP","CTP", ;S mean south,N mean North,EA mean East Asia
lats = (/ 22  , 30  /)
latn = (/ 30  , 45  /)
lonl = (/107  ,115  /)
lonr = (/130  ,130  /)
lev  = (/850,200/)
narea = dimsizes(area)
corr = new((/narea,nvar3,2*mxlag+1/),float)

do na = 0,narea-1,1
print(area(na))
;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(3)+".nc","r")
vars:= f->$var_name(3)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = 2
preci := new((/ncase,nyear,ntime,nlat,nlon/),float)
u     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
v     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
q     := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
q1    := new((/ncase,nyear,ntime,nlev,nlat,nlon/),float)
var   := new((/nvar2,ncase,nyear,ntime,nlat,nlon/),float)
corr1 := new((/nvar3,nyear,nlat,nlon,(mxlag+1)/),float)
corr2 := new((/nvar3,nyear,nlat,nlon,(mxlag+1)/),float)
corr11 := new((/nvar3,nyear,(mxlag+1)/),float)
corr22 := new((/nvar3,nyear,(mxlag+1)/),float)

if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
do nf = 0,1,1
do nt = 0,nyear-1,1
nv = 0
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
preci(nf,nt,:,:,:)   = f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
nv = 1
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
preci(nf,nt,:,:,:)   = preci(nf,nt,:,:,:) + f->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa

nv = 2
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})

nv = 3
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})

nv = 4
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q(nf,nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})

nv = 5
f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
q1(nf,nt,:,:,:,:) = tofloat(f->$var_name(nv)$(select_time,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)}))
end do
end do
u  = linmsg_n(u ,-1,5)
v  = linmsg_n(v ,-1,5)
q  = linmsg_n(q ,-1,5)
q1 = linmsg_n(q1,-1,5)
var(0,:,:,:,:,:) = preci*24*3600*1000 ;convert m/s to mm/day
delete([/preci/])
printVarSummary(u)

div := uv2dv_cfd( u, v, vars&lat, vars&lon, 2)
if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb    
    opt = False
    edd_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    edd_v = bw_bandpass_filter(v,fca,fcb,opt,1)
    edd_q = bw_bandpass_filter(q,fca,fcb,opt,1)
else
    edd_q := q - conform(q,dim_avg_n(q,1),(/0,2,3,4,5/))
    edd_u := u - conform(u,dim_avg_n(u,1),(/0,2,3,4,5/))
    edd_v := v - conform(v,dim_avg_n(v,1),(/0,2,3,4,5/))   
end if
delete([/u,v,q/])
EKE    := (edd_v*edd_v + edd_u*edd_u)/2.0   ;ps-200hPa integral get (nyear,ntime,nlat,nlon)
var(1,:,:,:,:,:) = uv2dv_cfd( edd_q(:,:,:,0,:,:)*edd_u(:,:,:,0,:,:), edd_q(:,:,:,0,:,:)*edd_v(:,:,:,0,:,:), vars&lat, vars&lon, 2) ;QTE850
delete([/edd_u,edd_v,edd_q/])
var(2,:,:,:,:,:) = EKE(:,:,:,0,:,:)  ;EKE850 :
var(3,:,:,:,:,:) = EKE(:,:,:,1,:,:)  ;EKE200 :
var(4,:,:,:,:,:) =  q1(:,:,:,0,:,:)  ;q1_850 :
var(5,:,:,:,:,:) =  q1(:,:,:,1,:,:)  ;q1_200 :
var(6,:,:,:,:,:) = div(:,:,:,0,:,:)  ;div850 :
var(7,:,:,:,:,:) = div(:,:,:,1,:,:)  ;div200 :
var(:,1,:,:,:,:) = var(:,1,:,:,:,:) - var(:,0,:,:,:,:)
delete([/q1/])

copy_VarMeta(vars(0,0,:,:),var(0,0,0,0,:,:))
;var := wgt_areaave(var,1.0,1.0,0) ;The rightmost dimensions should correspond to "latitude" (lat) and "longitude" (lon)
var!0 = "nvar"
var!1 = "case"
var!2 = "year"
var!3 = "time"

nv3 = 0
do nv1 = 0,nvar2-2, 1
do nv2 = nv1+1, nvar2-1, 1
    ;corr11(nv3,:,:) = esccr(var(nvar|nv1,case|fnc,year|:,time|:),var(nvar|nv2,case|fnc,year|:,time|:),mxlag)
    ;corr22(nv3,:,:) = esccr(var(nvar|nv2,case|fnc,year|:,time|:),var(nvar|nv1,case|fnc,year|:,time|:),mxlag)
    corr1(nv3,:,:,:,:) = esccr(var(nvar|nv1,case|fnc,year|:,lat|:,lon|:,time|:),var(nvar|nv2,case|fnc,year|:,lat|:,lon|:,time|:),mxlag)
    corr2(nv3,:,:,:,:) = esccr(var(nvar|nv2,case|fnc,year|:,lat|:,lon|:,time|:),var(nvar|nv1,case|fnc,year|:,lat|:,lon|:,time|:),mxlag)
    nv3 = nv3 + 1
end do
end do
print(nv3)

copy_VarMeta(var(0,0,:,:,:,:),corr1(:,:,:,:,0))
copy_VarMeta(var(0,0,:,:,:,:),corr2(:,:,:,:,0))
do nm = 0,mxlag,1
    corr11(:,:,nm) = wgt_areaave(corr1(:,:,:,:,nm),1.0,1.0,0) ;The rightmost dimensions should correspond to "latitude" (lat) and "longitude" (lon)
    corr22(:,:,nm) = wgt_areaave(corr2(:,:,:,:,nm),1.0,1.0,0)
end do
corr(na,:,0:mxlag-1) = dim_avg_n(corr22(:,:,1:mxlag:-1),1) ;calc climate and negative lag
corr(na,:,mxlag:)    = dim_avg_n(corr11(:,:,0:mxlag:),1)   ;calc climate and positive lag
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
fileatt@discription   = nvar3 + " delay correlation between preci,QTE850,EKE850,EKE200,q1_850,q1_200,div850,div200"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"area","nmber","delay"/);
dimSizes = (/narea,nvar3   ,2*mxlag+1/) ;
dimUnlim = (/False,False,False      /);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"area" ,"integer",(/"area" /))
filevardef(ncdf,"nmber","integer",(/"nmber"/))
filevardef(ncdf,"delay","integer",(/"delay"/))
filevardef(ncdf,"corr" ,"float"  ,(/"area","nmber","delay"/)) 

ncdf->narea  = ispan(1,narea,1)
ncdf->nmber  = ispan(1,nvar3,1)
ncdf->delay  = ispan(-mxlag,mxlag,1)
ncdf->corr   = (/corr/)
end

