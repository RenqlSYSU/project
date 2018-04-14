;***************************************************************************
;filename:      calc_preci_TP_daily_wavelet.ncl
;read the 27 years JJA daily PRECC and PRECL of two cases at East Asia(20N£­45N£¬107-130E)
;then do the wavelet analysis. 
;                                            by Ql_Ren
;                                           2018/03/30
;******************************************************************************
begin
area =(/"TP","eastAsia"/)
lats = (/25 , 20/) ;-20
latn = (/40 , 45/)  ;60
lonl = (/70 ,107/)  ;20
lonr = (/105,122/);130 ;220
na = 0;1

path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/wavelet-20180410/mdata/" + (/"CTRL","NUDG6h","NUDG24h"/) + "-preci_" + area(na) + "_JJAdaily_wavelet.nc"
case     = (/"CTRL","NUDG6h","NUDG24h"/)
var_name = (/"PRECC","PRECL"/)
years    = ispan(1979,2005,1)   ;the calculate time

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
files1     = systemfunc("ls " + filein(0) + "*" + var_name(0) + ".nc")
f          = addfiles(files1,"r")
ListSetType(f,"cat")
time        = f[:]->time
times       = cd_calendar(time,0)
select_time = ind(times(:,1).ge.6.and.times(:,1).le.8);.and.time(:,2).le.30)
ntime       = dimsizes(select_time)
print(ntime)

vars  := f[:]->$var_name(0)$(0,{lats(na):latn(na)},{lonl(na):lonr(na)})
nyear = dimsizes(years)
nvar  = dimsizes(var_name)
nfile = dimsizes(filein)
var_d = new((/ntime/),float)

do nf = 0, nfile-1,1
;--------------------------------------------------------------
;read the data
;--------------------------------------------------------------
var_d = 0
    do nv = 0, nvar-1,1 
        files1     = systemfunc("ls " + filein(nf) + "*" + var_name(nv) + ".nc")
        f          = addfiles(files1,"r")
        ListSetType(f,"cat")
        temp := f[:]->$var_name(nv)$(select_time,{lats(na):latn(na)},{lonl(na):lonr(na)})
        temp := wgt_areaave_Wrap(temp,1.0,1.0,0) 
        var_d = var_d + temp
    end do  ;the loop of PRECC and PRECL
var_d  = var_d*1000*3600*24  ;convert unit from m/s to mm/day

;************************************
; compute wavelet
;************************************
mother  = 0    ;0 = Morlet£¬1 = Paul£¬2 = DOG (derivative of Gaussian)
param   = 6.0  ;The mother wavelet parameter. If param < 0, then the default is used,
;Morlet k0 (wavenumber), default is 6. Paul m (order), default is 4. DOG m (m-th derivative), default is 2.

N       = ntime
dt      = 1    ;timesteps of time series, the unit here is day
s0      = dt   ;The smallest scale of the wavelet, which is typically is equal to 2*dt. s0 = dt for Morlet; s0 = dt/4 for Paul
dj      = 0.25 ;The spacing between discrete scales, which is typically equal to 0.25. A smaller value will give better scale resolution, but will be slower.
jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))  
;The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj].

npad    = N    ;The total number of points (including padding) to use for the wavelet transform
nadof   = 0    ;Currently ignored (set to zero).

noise   = 1    ;0 means use a white noise for significance test. 1 means use a red noise for significance test. (Most commonly, noise = 1.)
siglvl  = .05  ;The significance level to use
isigtest= 0    ; 0 means do a regular chi-square test, i.e. Eqn (18) from Torrence and Compo. 1 means do a "time-average" test on the global wavelet spectrum.

var = wavelet(var_d,mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
printVarSummary(var)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nf))
ncdf = addfile(fileout(nf),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the wavelet analysis of JJA daily rainfall in " + area(na)
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"part" ,"jtot","time"/)
dimSizes = (/  2    ,jtot  ,ntime /)
dimUnlim = (/False  ,False ,False /)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"part","integer"   ,(/"part"/))
filevardef(ncdf,"jtot","integer"   ,(/"jtot"/))
filevardef(ncdf,"time","integer"   ,(/"time"/))
filevardef(ncdf,"wave",typeof(var) ,(/"part","jtot","time"/))

;Define var Attribute
part    = (/0,1/)
part@long_nmae   = "the real (0,:,:) and imaginary parts (1,:,:) of the wavelet transform"
jtot@long_name   = "The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj]."
filevarattdef(ncdf,"part",part)
filevarattdef(ncdf,"jtot",jtot)
filevarattdef(ncdf,"time",time)
filevarattdef(ncdf,"wave",var)

ncdf->part   = (/0,1/)  
ncdf->jtot   = ispan(1,jtot,1) 
ncdf->time   = ispan(1,N,1)
ncdf->wave   = (/var/)
end do
end

