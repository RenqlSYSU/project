;*********************************************************************
;filename:      calc_area_avrg_ts.ncl
;output_file:   area_avrg_ts_winter.nc
;judge the ind of winer time(9,10,11,12,1,2) then declare an array ts(nvar,ntime)
;next read the Ruyuan rain and area 850 and 300 divg,850 vort, 700omega, flux, flux-divg
;                                            by Ql_Ren
;                                           2018/03/06
;*********************************************************************************
begin
path     = "/home/ys17-19/renql/project/2017downscaling"
var_name = (/"rain","omega","uwnd","vwnd","shum","pres.sfc"/)
filein1  = path + "/data/rain7016.nc"
filein2  = path + "/data/NECP_NCAR-daily/" + var_name(1:5) + "/" + var_name(1:5) + "*"
fileout  = path + "/regression/mdata/" + "area_avrg_ts_winter.nc" 
fileout2 = path + "/regression/mdata/" + "area_avrg_ts_winter.dat" 
fileout3 = path + "/regression/mdata/" + "area_avrg_ts_winter.txt" 
years    = ispan(1970,2016,1)   ;the calculate time
season   = "Winter"
g        = 9.8 ;m/(s*s)

;area of field 
ts_var   = (/"rain","700omega","850divg","300divg","850vort","flux","flux-divg"/)
lats     = (/   0  ,    22    ,    22   ,    25   ,    20   ,  22  ,     20    /)
latn     = (/   0  ,    27    ,    27   ,    28   ,    25   ,  27  ,     25    /)
lonl     = (/   0  ,    110   ,    110  ,    112  ,    107  ,  110 ,     108   /)
lonr     = (/   0  ,    117   ,    117  ,    120  ,    111  ,  117 ,     113   /)
lev      = (/   0  ,    700   ,    850  ,    300  ,    850  /)

files       = systemfunc("ls " + filein2(0)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
time        = f2[:]->time
time_cd     = cd_calendar(time,0)
select_time = ind(time_cd(:,1).ge.9.or.time_cd(:,1).le.2)     ;winter
;select_time = ind(time_cd(:,1).ge.3.and.time_cd(:,1).le.8)   ;summer
;select_time = ind(time_cd(:,1).ge.1.and.time_cd(:,1).le.12)  ;all year
ntime       = dimsizes(select_time)
nvar        = dimsizes(ts_var)
ts          = new((/nvar,ntime/),float)
print(ntime)
printVarSummary(select_time)

;==========================================================================
;read the winter(9,10,11,12,1,2) data,2016-1948=68,1970-1948=22 and ruyuan rainfall
;===========================================================================
f1      = addfile(filein1,"r")
ts(0,:) = f1->$var_name(0)$(select_time)  ;get var1(time)
print("read "+ts_var(0))

;----700omega---------
nv = 1
files       = systemfunc("ls " + filein2(nv-1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var        := f2[:]->$var_name(nv)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})
ts(nv,:)    = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv+" lev = "+lev(nv))
printVarSummary(var)

;----850divg---------
nv = 2
files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_u      := f2[:]->$var_name(2)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v      := f2[:]->$var_name(3)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

var        := uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 2)
copy_VarMeta(var_u,var)
ts(nv,:)    = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv+" lev = "+lev(nv))
printVarSummary(var_u)
printVarSummary(var_v)

;------300divg---------
nv = 3
files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_u      := f2[:]->$var_name(2)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v      := f2[:]->$var_name(3)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

var        := uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 2)
copy_VarMeta(var_u,var)
ts(nv,:)    = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv+" lev = "+lev(nv))
printVarSummary(var_u)
printVarSummary(var_v)

;------850vort-----------
nv = 4
files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_u      := f2[:]->$var_name(2)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v      := f2[:]->$var_name(3)$(select_time,{lev(nv)},{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

var        := uv2vr_cfd(var_u, var_v, var_v&lat, var_v&lon, 2)
copy_VarMeta(var_u,var)
ts(nv,:)    = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv+" lev = "+lev(nv))
printVarSummary(var_u)
printVarSummary(var_v)

;----------flux----------
nv = 5
files       = systemfunc("ls " + filein2(3)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_shum   := f2[:]->$var_name(4)$(select_time,:,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})
var_shum    = var_shum*1000  ;convert unit from kg/kg to g/kg

files       = systemfunc("ls " + filein2(4)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_pres   := f2[:]->pres(select_time,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})
var_pres    = var_pres/100   ;convert unit from Pa to hPa

var        := wgt_vert_avg_beta ( var_shum&level, var_shum, var_pres, 0, 0 ) / g 
;the first 0 means the units of p and psfc are hPa, the second 0 means the pressure-thickness weighted vertical sum
copy_VarMeta(var_pres,var)
ts(nv,:)    = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv)
printVarSummary(var_shum)
printVarSummary(var_pres)

;-------flux-divg-----------
nv = 6
files       = systemfunc("ls " + filein2(3)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_shum   := f2[:]->$var_name(4)$(select_time,:,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})
var_shum    = var_shum*1000  ;convert unit from kg/kg to g/kg

files       = systemfunc("ls " + filein2(4)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_pres   := f2[:]->pres(select_time,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})
var_pres    = var_pres/100   ;convert unit from Pa to hPa

files       = systemfunc("ls " + filein2(1)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_u      := f2[:]->$var_name(2)$(select_time,0:7,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

files       = systemfunc("ls " + filein2(2)) 
f2          = addfiles(files(22:68),"r")
ListSetType(f2,"cat")
var_v      := f2[:]->$var_name(3)$(select_time,0:7,{lats(nv):latn(nv)},{lonl(nv):lonr(nv)})

var_u    = var_u*var_shum/g
var_v    = var_v*var_shum/g
var_divg = uv2dv_cfd(var_u, var_v, var_v&lat, var_v&lon, 3) ;calculate the divergence
copy_VarMeta(var_u,var_divg)
var     := wgt_vert_avg_beta ( var_shum&level, var_divg, var_pres, 0, 0 )
ts(nv,:) = wgt_areaave(var,1,1,0)
print("read "+ts_var(nv)+",nv = "+nv)
printVarSummary(var_shum)
printVarSummary(var_pres)

;========================================================================
;save the data
;=======================================================================
fbindirwrite(fileout2,ts)

;write_table(fileout3,"w",ts_var,"%s%s%s%s%s%s%s")
time_YYYYMMDD = cd_calendar(time(select_time),-2)
alist = (/time_YYYYMMDD,ts(0,:),ts(1,:),ts(2,:),ts(3,:),ts(4,:),ts(5,:),ts(6,:)/)
write_table(fileout3,"w",ts,"%8d%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f")

system("rm -f " + fileout)
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribute
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "time series of rain,700omega,850divg,300divg,850vort,flux,flux-divg" + \ 
                        " in winter(9,10,11,12,1,2) from 1970 to 2016"
                        ;" from 1970 to 2016"
                        ;" in summer(3,4,5,6,7,8) from 1970 to 2016"
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"var","time"/)
dimSizes = (/nvar , ntime/)
dimUnlim = (/False, False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"var" ,"integer"   ,(/"var"/))
filevardef(ncdf,"time",typeof(time),getvardims(time))
filevardef(ncdf,"ts"  ,typeof(ts)  ,(/"var","time"/))

;Define var Attribute
filevarattdef(ncdf,"time",time)
ts@long_name = "the units of seven variables is mm, Pa/s, /s, /s, /s, g/(m*m), g/(s*m*m)"
filevarattdef(ncdf,"ts",ts)

ncdf->var  = (/0,1,2,3,4,5,6/) 
ncdf->time = (/time(select_time)/)
ncdf->ts      = (/ts/)
end
