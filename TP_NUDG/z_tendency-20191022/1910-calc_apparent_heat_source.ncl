;***************************************************************************
;use the post-processing data to calculate daily apparent heat source
;1. do a loop of year and read U, V, T, W
;2. calc Q1 and store it
;
;The ncl functions have been used are: 
;   center_finite_diff_n
;                                            by Ql_Ren
;                                           2019/10/23
;******************************************************************************
begin
path    = "/home/ys17-19/renql"
case_option = 2 ;0 means NUDG experiment; 1 means noTP experiment

if(case_option.eq.0) then
    case    = (/"CTRL","NUDG24h","NUDG6h"/)
    filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
    years   = ispan(1979,2005,1)   ;the calculate time
    fileout = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                      "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                      "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/); + "Q1.nc"
else
    ;remove TP experiment, sample number is ten years 
    case    = (/"TP_CTRL","TP_CR"/) ;TP_CTRL is 0021-0030; TP_CR is 0041-0050
    filein  = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
    years   = ispan(21,30,1)   ;the calculate time
    fileout = path + "/model/" + case + "/pro/" + case + ".cam.h1.YEAR.00" ; + year + ".daily." + var_name + ".nc"
end if

var_name= (/"U","V","T","OMEGA"/)
ncase   = dimsizes(case)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats = -90 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;
lev  = (/1000,925,850,800,750,700,650,600,550,500,450,400,350,300,250,200,150,100/) ;hPa
month_s = 1 
month_e = 12

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein(0)+years(0)+".daily."+var_name(0)+".nc","r")
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
u     = new((/ntime,nlev,nlat,nlon/),float)

iday  = 0
ilev  = 1
ilat  = 2
ilon  = 3
lat  = f->lat
lon  = f->lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
coslat = conform(u,cos(lat*pi/180.0),ilat)
lev4d  = conform(u, lev*100, ilev)
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
loncyc = True  ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nf = 1,1,1
if((case_option.eq.1).and.(nf.eq.1))
    years   = ispan(41,50,1)   ;the calculate time
end if

do nt = 0,nyear-1,1
print("handle with " + case(nf) + " " + years(nt))
;==========================================================
;read the data and calc apparent heat source Q1
;============================================================
nv = 0
f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
u  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 1
f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
v  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 2
f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
t  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})

nv = 3
f  = addfile(filein(nf) + years(nt) + ".daily." + var_name(nv) + ".nc","r")
w  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
printVarSummary(w)

sigm0 = (R/cp)*(t/lev4d) - center_finite_diff_n(t,lev*100,cyclic,opt,ilev)
Q1    =   (u/a/coslat) * center_finite_diff_n(t,lon,loncyc,opt,ilon) \
        + (v/a       ) * center_finite_diff_n(t,lat,cyclic,opt,ilat) \
        - sigm0 * w + center_finite_diff_n(t,24*3600,cyclic,opt,iday)
delete([/v,t,w/])
Q1    = Q1*cp
print(Q1(15,5,30,30))

Q1@long_name = "Apparent Heat Source calc from temperature equation"
Q1@units     = "W/kg" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+years(nt)+".daily.Q1.nc")
ncdf = addfile(fileout(nf)+years(nt)+".daily.Q1.nc","c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "Daily apparent heat source (W/kg) calculated from the temperature equation" 
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
filevardef(ncdf,"Q1"  ,typeof(Q1  ),(/"time","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",u&time)
filevarattdef(ncdf,"lev",u&level)
filevarattdef(ncdf,"lat",u&lat)
filevarattdef(ncdf,"lon",u&lon)
filevarattdef(ncdf,"Q1" ,Q1   )

ncdf->time = (/u&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/u&lat/)
ncdf->lon  = (/u&lon/)
ncdf->Q1   = (/Q1   /)
end do
end do
end

