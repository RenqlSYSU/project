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
filein  = "/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_Q1_1979-2020.nc"
fileout = "/home/ys17-23/Extension2/renql/ERA5_mon/ERA5_mon_intQ1_1979-2020_sm9.nc"

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
nlev  = dimsizes(lev)

;==========================================================
;read the data and calc apparent heat source Q1
;============================================================
f   = addfile(filein,"r")
vars= f->Q1
printVarSummary(vars)
ntime = dimsizes(vars&time)
nlat  = dimsizes(vars&latitude)
nlon  = dimsizes(vars&longitude)

f  = addfile("/home/ys17-23/Extension/ERA5/ERA5-monthly/surface/sur_pressure_mon_1x1_1979_2020.nc","r")
ps = short2flt(f->sp(:,{lats:latn},{lonl:lonr}))
dp = dpres_plevel_Wrap(lev*100,ps,10000,1)  ;unit is Pa, get dp (year,time,lev,lat,lon)
Q1 = dim_sum_n(vars*dp,1)/g
delete([/ps,dp/])
Q1 = smth9(Q1,0.5,0.25,True)
Q1@long_name = "Vertically integrated Apparent Heat Source calc from temperature equation"
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
fileatt@discription   = "monthly apparent heat source (W/kg) calculated from the temperature equation" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"time","latitude","longitude"/);
dimSizes = (/ntime ,nlat ,nlon/) ;
dimUnlim = (/False ,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"time",typeof(vars&time),getvardims(vars&time))
filevardef(ncdf,"latitude" ,typeof(vars&latitude),(/"latitude"/))
filevardef(ncdf,"longitude" ,typeof(vars&longitude),(/"longitude"/))
filevardef(ncdf,"Q1"  ,typeof(Q1  ),(/"time","latitude","longitude"/))

;define var attribute
filevarattdef(ncdf,"time",vars&time)
filevarattdef(ncdf,"latitude",vars&latitude)
filevarattdef(ncdf,"longitude",vars&longitude)
filevarattdef(ncdf,"Q1" ,Q1   )

ncdf->time = (/vars&time/)
ncdf->latitude  = (/vars&latitude/)
ncdf->longitude  = (/vars&longitude/)
ncdf->Q1   = (/Q1   /)
end

