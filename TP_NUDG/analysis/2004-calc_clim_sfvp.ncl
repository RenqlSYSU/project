;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 0 ;0 means NUDG experiment; 1 means noTP experiment
nmonth = 4
month_s = 1
month_e = 8

if(case_option.eq.0) then
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + (/"CTRL","NUDG24h","NUDG6h"/) + "-Clim_2month_ave_sfvp.nc"
    years    = ispan(1979,2005,1)   ;the calculate time
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    fileout = path + "/project/TP_NUDG/analysis/mdata/" + case + "-Clim_2month_ave_sfvp.nc"
    years   = ispan(21,30,1)   ;the calculate time
end if
var_name = (/"U","V"/)

;area: east Asia
lats = -90 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,800,750, 700,650,600,550,500, 450,400,350,300,250, 200,150,100/) ;18 levels ,hPa U

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
vars:= f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
time  = cd_calendar(f->time,0)
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
ncase = dimsizes(filein)

;==========================================================
;read the data
;============================================================
do nf = 0 ,2, 2
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if

u_m  = new((/nmonth,nlev,nlat,nlon/),float)
v_m  = new((/nmonth,nlev,nlat,nlon/),float)
do nm = 0, nmonth-1, 1
month_s = 2*nm + 1
month_e = 2*nm + 2
select_time := ind(time(:,1).ge.month_s.and.time(:,1).le.month_e)
ntime = dimsizes(select_time)
print("===================================================================")
print("month from " + month_s + " to " + month_e + ", ntime = " + ntime)
u  := new((/nyear,ntime,nlev,nlat,nlon/),float)
v  := new((/nyear,ntime,nlev,nlat,nlon/),float)

do nt = 0,nyear-1,1
    print("read: " + years(nt))
    nv = 0
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    u(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f  = addfile(filein(nf)+years(nt)+".daily."+var_name(nv)+".nc","r")
    v(nt,:,:,:,:) = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end do
u := dim_avg_n_Wrap(u,0)
v := dim_avg_n_Wrap(v,0)
u_m(nm,:,:,:) = dim_avg_n_Wrap(u,0) 
v_m(nm,:,:,:) = dim_avg_n_Wrap(v,0) 
end do
delete([/u,v/])
u_m = linmsg_n(u_m,-1,3)
v_m = linmsg_n(v_m,-1,3)
u_m = linmsg_n(u_m,-1,2)
v_m = linmsg_n(v_m,-1,2)
u_m = linmsg_n(u_m,-1,1)
v_m = linmsg_n(v_m,-1,1)

if(any(ismissing(u_m))) then
    print("uwnd has missing")
end if
if(any(ismissing(v_m))) then
    print("vwnd has missing")
end if

vort = new((/nmonth,nlev,nlat,nlon/),float)
divg = new((/nmonth,nlev,nlat,nlon/),float)
sf   = new((/nmonth,nlev,nlat,nlon/),float)
vp   = new((/nmonth,nlev,nlat,nlon/),float)
uv2vrdvf(u_m,v_m,vort,divg) ;Computes the vorticity and divergence via spherical harmonics
uv2sfvpf(u_m,v_m,sf,vp) ;calculate stream function and velocity potential via spherical harmonics
delete([/u_m,v_m/])
ur   = new((/nmonth,nlev,nlat,nlon/),float)
vr   = new((/nmonth,nlev,nlat,nlon/),float)
ud   = new((/nmonth,nlev,nlat,nlon/),float)
vd   = new((/nmonth,nlev,nlat,nlon/),float)
vr2uvf(vort,ur,vr)    ;calculate the rotational wind components via spherical harmonics
dv2uvf(divg,ud,vd)    ;calculate the divergent wind components via spherical harmonics

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "all variables are calculated by spherical harmonics" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"month","lev","lat","lon"/);
dimSizes = (/nmonth ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vort",typeof(vort),(/"month","lev","lat","lon"/))
filevardef(ncdf,"divg",typeof(divg),(/"month","lev","lat","lon"/))
filevardef(ncdf,"sf",typeof(sf),(/"month","lev","lat","lon"/))
filevardef(ncdf,"vp",typeof(vp),(/"month","lev","lat","lon"/))
filevardef(ncdf,"ur",typeof(ur),(/"month","lev","lat","lon"/))
filevardef(ncdf,"vr",typeof(vr),(/"month","lev","lat","lon"/))
filevardef(ncdf,"ud",typeof(ud),(/"month","lev","lat","lon"/))
filevardef(ncdf,"vd",typeof(vd),(/"month","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)

ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->vort = (/vort/)
ncdf->divg = (/divg/)
ncdf->sf = (/sf/)
ncdf->vp = (/vp/)
ncdf->ur = (/ur/)
ncdf->vr = (/vr/)
ncdf->ud = (/ud/)
ncdf->vd = (/vd/)
end do
end

