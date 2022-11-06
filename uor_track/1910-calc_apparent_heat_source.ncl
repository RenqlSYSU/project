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
filein  = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/"
fileout = "/gws/nopw/j04/ncas_generic/users/renql/ERA5_subdaily/Q1/ERA5_";Q1_2012.nc
years = ispan(1980,1989,1)

var_name= (/"U","V","T","OMEGA"/)
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
;lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200,150,100/) ;hPa
month_s = 1 
month_e = 12

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein+"u/ERA5_u_1989.nc","r")
vars:= f->$var_name(0)$
time = cd_calendar(f->time,0)
ntime = dimsizes(time)
nlat  = dimsizes(vars&latitude)
nlon  = dimsizes(vars&longitude)
nlev  = dimsizes(vars&level)
nyear = dimsizes(years)
nvar  = dimsizes(var_name)

iday  = 0
ilev  = 1
ilat  = 2
ilon  = 3
lev  = vars&level
lat  = vars&latitude
lon  = vars&longitude
lat  = lat*pi/180.0
lon  = lon*pi/180.0
coslat = conform(vars, cos(lat),ilat)
lev4d  = conform(vars, lev*100, ilev)
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
loncyc = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

do nt = 0,nyear-1,1
print("handle with " + years(nt))
;==========================================================
;read the data and calc apparent heat source Q1
;============================================================
nv = 0
f  = addfile(filein+var_name(nv)+"/ERA5_"+var_name(nv)+"_"+years(nt)+".nc","r")
u  = f->$var_name(nv)$
nv = 1
f  = addfile(filein+var_name(nv)+"/ERA5_"+var_name(nv)+"_"+years(nt)+".nc","r")
v  = f->$var_name(nv)$
nv = 2
f  = addfile(filein+var_name(nv)+"/ERA5_"+var_name(nv)+"_"+years(nt)+".nc","r")
t  = f->$var_name(nv)$
nv = 3
f  = addfile(filein+var_name(nv)+"/ERA5_"+var_name(nv)+"_"+years(nt)+".nc","r")
w  = f->$var_name(nv)$

files := systemfunc("ls "+filein(nf)+"*.daily.PS.nc" )
f      = addfile(files(nt),"r")
ps   = f->PS(select_time,{lats:latn},{lonl:lonr}) ;get an array (ncase,nmonth,nlat,nlon)
psn  = conform(u,ps,(/0,2,3/))
delete(ps)
levn = conform(u,lev*100,(/1/))
u  = mask( u, psn.ge.levn,True)
v  = mask( v, psn.ge.levn,True)
t  = mask( t, psn.ge.levn,True)
w  = mask( w, psn.ge.levn,True)
delete([/psn,levn/])

sigm0 = (R/cp)*(t/lev4d) - center_finite_diff_n(t,lev*100,cyclic,opt,ilev)
Q1    =   (u/a/coslat) * center_finite_diff_n(t,lon,loncyc,opt,ilon) \
        + (v/a       ) * center_finite_diff_n(t,lat,cyclic,opt,ilat) \
        - sigm0 * w + center_finite_diff_n(t,24*3600,cyclic,opt,iday)
delete([/u,v,t,w,sigm0/])
Q1    = Q1*cp
print(Q1(15,5,30,30))

Q1@long_name = "Apparent Heat Source calc from temperature equation"
Q1@units     = "W/kg" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf)+years(nt)+".daily.Q11.nc")
ncdf = addfile(fileout(nf)+years(nt)+".daily.Q11.nc","c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "Daily apparent heat source (W/kg) calculated from the temperature equation, mask value bellow topography" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","lev","lat","lon"/);
dimSizes = (/ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"Q1"  ,typeof(Q1  ),(/"time","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"Q1" ,Q1   )

ncdf->time = (/vars&time/)
ncdf->lev  = (/lev/)
ncdf->lat  = (/vars&lat/)
ncdf->lon  = (/vars&lon/)
ncdf->Q1   = (/Q1   /)
end do
end do
end

