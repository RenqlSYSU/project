;***************************************************************************
;use the post-processing data to calculate climatology daily transient eddies
;1. remove annual cycle;  2. calc transient eddies kinetic energy (EkE)
;3. calc climatology daily EKE
;The ncl functions have been used are: dpres_plevel_Wrap,dim_sum_n_Wrap,dim_avg_n_Wrap
;dim_rmvmean_n_Wrap
;                                            by Ql_Ren
;                                           2019/01/25
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment
nmonth = 4
month_s = 1
month_e = 8

if(case_option.eq.0) then
    case    = (/"CTRL","NG6h","NG24h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-Clim_2month_ave_advT.nc"
end if  
if(case_option.eq.1) then
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_2month_ave_advT.nc"
end if
if(case_option.eq.2) then
    case    = (/"F2000","NGTR"/) ;
    filein  = path +(/"/model/F2000_CAM5/F2000_CAM5.cam.h1.YEAR.00",\
                      "/model/F2000_CAM5_TP_NUDG/F2000_CAM5_TP_NUDG.cam.h1.YEAR.00"/)
    years   = ispan(1,30,1)   ;the calculate time
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_2month_ave_advT.nc"
end if

var_name= (/"U","V","T"/)
ncase   = dimsizes(case)

;area: east Asia
lats =  0  ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/850,500,200/) 
;lev  = (/1000,925,850,700,600,  500,400,350,300,250,  200,150,100/) 

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(9)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
time = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)
var_month  = new((/3,nyear,nmonth,nlev,nlat,nlon/),double)
month_ave  = new((/3,      nmonth,nlev,nlat,nlon/),double)
month_vari = new((/3,      nmonth,nlev,nlat,nlon/),double) ;be used to exam the significance of the mean difference between the two groups

do nf = 0,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if
;==========================================================
;read the data
;============================================================
u  := new((/nyear,ntime,nlev,nlat,nlon/),float)
v  := new((/nyear,ntime,nlev,nlat,nlon/),float)
t  := new((/nyear,ntime,nlev,nlat,nlon/),float)
do nt = 0,nyear-1,1
nv = 0
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
files := systemfunc("ls "+filein(nf)+"*.daily."+var_name(nv)+".nc" )
f  = addfile(files(nt),"r")
print(" handle with " + files(nt))
;f      = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
printVarSummary(u)

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)
iyear = 0
itime = 1
ilev  = 2
ilat  = 3
ilon  = 4
lat  = vars&lat
lon  = vars&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
opt    = 0     ;used by center_finite_diff_n, no meanging 
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 
coslat5d = conform(u,cos(lat),ilat)
var = -(u/a/coslat5d)*center_finite_diff_n(t,lon,cyclic,opt,ilon) \
      -(v/a         )*center_finite_diff_n(t,lat,cyclic,opt,ilat)
var@long_name = "temperature advection: total, stationary, transient"
var@units     = "K/s"
delete(coslat5d)

time1 = cd_calendar(f->time(select_time),0)
do nm = 0, nmonth-1, 1
    mons = 2*nm + 1
    mone = 2*nm + 2
    print("month from " + mons + " to " + mone) 
    select_time1 := ind(time1(:,1).ge.mons.and.time1(:,1).le.mone)
    var_month(0,:,nm,:,:,:) = dim_avg_n_Wrap(var(:,select_time1,:,:,:),itime) ;get (nyear,nmonth,nlat,nlon)
    u_month = dim_avg_n_Wrap(u(:,select_time1,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    v_month = dim_avg_n_Wrap(v(:,select_time1,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    t_month = dim_avg_n_Wrap(t(:,select_time1,:,:,:),1) ;get (nyear,nmonth,nlat,nlon)
    coslat4d = conform(u_month,cos(lat),ilat-1)
    var_month(1,:,nm,:,:,:) = -(u_month/a/coslat4d)*center_finite_diff_n(t_month,lon,cyclic,opt,ilon-1) \
                              -(v_month/a         )*center_finite_diff_n(t_month,lat,cyclic,opt,ilat-1)
    delete(coslat4d)
    u_edd = u(:,select_time1,:,:,:)-conform(u(:,select_time1,:,:,:),u_month,(/0,2,3,4/))
    v_edd = v(:,select_time1,:,:,:)-conform(v(:,select_time1,:,:,:),v_month,(/0,2,3,4/))
    t_edd = t(:,select_time1,:,:,:)-conform(t(:,select_time1,:,:,:),t_month,(/0,2,3,4/))
    delete([/u_month,v_month,t_month/])
    coslat5d = conform(u_edd,cos(lat),ilat)
    term  = -(u_edd/a/coslat5d)*center_finite_diff_n(t_edd,lon,cyclic,opt,ilon) \
            -(v_edd/a         )*center_finite_diff_n(t_edd,lat,cyclic,opt,ilat)
    var_month(2,:,nm,:,:,:) = dim_avg_n_Wrap(term,1)
    delete([/u_edd,v_edd,t_edd,term,coslat5d/])
end do
delete([/u,v,t/])
month_ave  = dim_avg_n_Wrap(var_month,1) ;get (nmonth,nlat,nlon)
month_vari = dim_variance_n_Wrap(var_month,1) ;get (nmonth,nlat,nlon)
printVarSummary(month_ave)
print(month_ave(0,1,2,10,10:20))

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "rmvannual cycle"
fileattdef(ncdf,fileatt)  

;define all coordinate`
dimNames = (/"nvar","month","lev","lat","lon"/);
dimSizes = (/   3  ,nmonth ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"month_ave"  ,typeof(month_ave) ,(/"nvar","month" ,"lev","lat","lon"/))
filevardef(ncdf,"month_vari" ,typeof(month_vari),(/"nvar","month" ,"lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"month_ave",var)
filevarattdef(ncdf,"month_vari",var)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->month_ave   = (/month_ave/)
ncdf->month_vari  = (/month_vari/)
end do
end

