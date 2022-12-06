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
filein  = "/home/ys17-23/Extension/ERA5/ERA5-monthly/pressure/"
fileout = "/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_Q1_1979-2020.nc"

var_name= (/"u","v","t","w"/)
g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0
a  = 6378388 ;the radius of earth, m
pi = atan(1.0)*4
print(pi)

;area: east Asia
lats =  0 ;
latn =  90 ;
lonl =  0  ;
lonr = 360 ;
lev  = (/1000,925,850,700,600,  500,400,350,300,250, 200,150,100/) ;hPa

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f    = addfile(filein+"u_mon_1x1_1979_2020.nc","r")
vars:= short2flt(f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr}))
ntime = dimsizes(vars&time)
nlat  = dimsizes(vars&latitude)
nlon  = dimsizes(vars&longitude)
nlev  = dimsizes(lev)
nvar  = dimsizes(var_name)

iday  = 0
ilev  = 1
ilat  = 2
ilon  = 3
lat  = vars&latitude
lon  = vars&longitude
lat  = lat*pi/180.0
lon  = lon*pi/180.0
coslat = conform(vars, cos(lat),ilat)
lev4d  = conform(vars, lev*100, ilev)
opt    = 0     ;used by center_finite_diff_n, no meanging 
cyclic = False ;used by center_finite_diff_n
loncyc = True ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

;==========================================================
;read the data and calc apparent heat source Q1
;============================================================
nv = 0
f  = addfile(filein+var_name(nv)+"_mon_1x1_1979_2020.nc","r")
u  = short2flt(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
nv = 1
f  = addfile(filein+var_name(nv)+"_mon_1x1_1979_2020.nc","r")
v  = short2flt(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
nv = 2
f  = addfile(filein+var_name(nv)+"_mon_1x1_1979_2020.nc","r")
t  = short2flt(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))
nv = 3
f  = addfile(filein+var_name(nv)+"_mon_1x1_1979_2020.nc","r")
w  = short2flt(f->$var_name(nv)$(:,{lev},{lats:latn},{lonl:lonr}))

;f  = addfile("/home/ys17-23/Extension/ERA5/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2020.nc","r")
;ps = short2flt(f->sp(:,{lats:latn},{lonl:lonr}))
;psn  = conform(u,ps,(/0,2,3/))
;delete(ps)
;levn = conform(u,lev*100,(/1/))
;u  = mask( u, psn.ge.levn,True)
;v  = mask( v, psn.ge.levn,True)
;t  = mask( t, psn.ge.levn,True)
;w  = mask( w, psn.ge.levn,True)
;delete([/psn,levn/])

sigm0 = (R/cp)*(t/lev4d) - center_finite_diff_n(t,lev*100,cyclic,opt,ilev)
Q1    =   (u/a/coslat) * center_finite_diff_n(t,lon,loncyc,opt,ilon) \
        + (v/a       ) * center_finite_diff_n(t,lat,cyclic,opt,ilat) \
        - sigm0 * w ;+ center_finite_diff_n(t,24*3600,cyclic,opt,iday)
delete([/u,v,t,w,sigm0/])
Q1    = Q1*cp
print(Q1(15,5,30,30))

Q1@long_name = "Apparent Heat Source calc from temperature equation"
Q1@units     = "W/kg" 

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "monthly apparent heat source (W/kg) calculated from the temperature equation, mask value bellow topography" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","level","latitude","longitude"/);
dimSizes = (/ntime ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"level" ,typeof(vars&level),(/"level"/))
filevardef(ncdf,"latitude" ,typeof(vars&latitude),(/"latitude"/))
filevardef(ncdf,"longitude" ,typeof(vars&longitude),(/"longitude"/))
filevardef(ncdf,"Q1"  ,typeof(Q1  ),(/"time","level","latitude","longitude"/))

;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"level",vars&level)
filevarattdef(ncdf,"latitude",vars&latitude)
filevarattdef(ncdf,"longitude",vars&longitude)
filevarattdef(ncdf,"Q1" ,Q1   )

ncdf->time = (/vars&time/)
ncdf->level  = (/vars&level/)
ncdf->latitude  = (/vars&latitude/)
ncdf->longitude  = (/vars&longitude/)
ncdf->Q1   = (/Q1   /)
end

