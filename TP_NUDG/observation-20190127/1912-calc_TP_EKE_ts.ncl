;******************************************************************************
;1. calculate TP area average vertically integrated EKE
;2. calculate the correlation between
;                                            by Ql_Ren
;                                           2019/08/19
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/TP_NUDG/observation-20190127/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
path4   = "/home/ys17-19/data-observation/ERA-interim/surface/pres/sp.interim.nc"
gtopfile = "/home/ys17-19/data-observation/ERA-interim/gtopo_1.5x1.5.nc" ;in order to select TP area
var_name = (/"uwnd","vwnd"/)
option   = 0 ;0 means remove annual cycle, 1 means bandpass filter
na = 0

;parameter for bw_bandpass_filter
ca = 10.0 ;day
cb = 2.5 

g  = 9.8 ;m/(s*s)
years    = ispan(1979,2016,1)
years@calendar = "standard"

;area of index, India,Tibet,Philippines,Maritime_continent
area = (/"TP"/)
lats = (/25 /)
latn = (/40 /)
lonl = (/70 /)
lonr = (/105/)
lev  = (/1000,925,850,700,600, 500,400,350,300, 250,200,175,150, 125,100/)
;lev  = (/1000,975,950,925,900, 875,850,825,800,775, 750,700,650,600,550, 500,450,400,350,300, 250,225,200,175,150, 125,100/)

if(option.eq.0) then
    fileout  = path2 + "mdata/ERA-Interim_intEKE_" + area(na) + "_ts.nc"
else
    fileout  = path2 + "mdata/ERA-Interim_intEKE_" + area(na) + "_ts-" + cb + "-" + ca + "day.nc"
end if

;==================================================================================
;read data
;=============================================================================
f_model = addfile(path3,"r")
time = f_model->time
printVarSummary(time)

f = addfile(gtopfile,"r")
phis   = short2flt(f->z(0,{lats(na):latn(na)},{lonl(na):lonr(na)}))
phis   = phis/9.8   ;convert unit from m2/s2 to gpm

f     = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc","r")
vars  = f->$var_name(0)$(0,0,{lats(na):latn(na)},{lonl(na):lonr(na)})
ntime = 365 
nyear = dimsizes(years)
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
nlev = dimsizes(lev)
u    = new((/nyear,ntime,nlev,nlat,nlon/),float)
v    = new((/nyear,ntime,nlev,nlat,nlon/),float)
dp   = new((/nyear,ntime,nlev,nlat,nlon/),float)
ps   = new((/nyear,ntime,nlat,nlon/),float)
ts   = new((/nyear,ntime/),float)
delete(vars)

;==========================================================
;calculate the EKE index and read the field variable
;============================================================
index  = new((/ntime/),integer)
index1 = new((/ntime/),integer)
index(0:58)   = ispan(0,58,1)
index(59:364) = ispan(60,365,1)

f_sfcp  = addfile(path4,"r")
start = 0
do nt = 0,nyear-1,1
print("read: " + years(nt) + ", start = " + start)
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(index,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    space = 366
    index1(0:58)   = ispan(start,(start+58),1)
    index1(59:364) = ispan((start+60),(start+space-1),1)
    ps(nt,:,:,:)   = f_sfcp->sp(index1,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(:,{lev},{lats(na):latn(na)},{lonl(na):lonr(na)})
    space = 365
    ps(nt,:,:,:)   = f_sfcp->sp(start:(start+space-1):1,{lats(na):latn(na)},{lonl(na):lonr(na)})   ;ps(time,lat,lon),double, surface pressure, unit is Pa
end if
start = start + space
dp(nt,:,:,:,:) = dpres_plevel_Wrap(lev*100,ps(nt,:,:,:),10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
end do
printVarSummary(ps)
printVarSummary(u)
printVarSummary(dp)
delete(ps)

if(option.eq.1) then
    fca = 1.0/ca
    fcb = 1.0/cb
    opt = False
    filt_u = bw_bandpass_filter(u,fca,fcb,opt,1)
    filt_v = bw_bandpass_filter(v,fca,fcb,opt,1)
    var   = dim_sum_n((filt_v*filt_v + filt_u*filt_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
else
    edd_u = u - conform(u,dim_avg_n(u,0),(/1,2,3,4/))
    edd_v = v - conform(v,dim_avg_n(v,0),(/1,2,3,4/))
    delete([/u,v/])
    var   = dim_sum_n((edd_v*edd_v + edd_u*edd_u)*dp,2)/g/2.0   ;get (nyear,ntime,nlat,nlon)
    delete([/edd_u,edd_v,dp/])
end if
var = mask(var,phis.gt.1500,True) 
ts  = wgt_areaave_Wrap(var,1.0,1.0,0)
ts@long_name = "Climatology daily TP area average vertically integrated EKE"
ts@units     = "J" 
printVarSummary(ts)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the TP area average EKE ts" 
fileattdef(ncdf,fileAtt)  

;define all coordinate
dimNames = (/"year","time"/);
dimSizes = (/nyear , 365  /) ;
dimUnlim = (/False ,False /);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"year" ,"integer",(/"year" /))
filevardef(ncdf,"time" ,typeof(time),(/"time"/))
filevardef(ncdf,"EKEts",typeof(ts),(/"year","time"/))

;define var attribute
filevarattdef(ncdf,"time" ,time)
filevarattdef(ncdf,"EKEts",ts)

ncdf->year = years
ncdf->time = time
ncdf->EKEts = (/ts/)
end

