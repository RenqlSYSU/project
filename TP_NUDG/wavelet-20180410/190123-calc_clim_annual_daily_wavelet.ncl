;***************************************************************************
;read var(nyear,nday),remove annual cycle and standardization
;then calc the wavelet for every year
;then average power,phase,signif,gws over all years.
;                                            by Ql_Ren
;                                           2018/03/30
;******************************************************************************
begin
area =(/"STP","CTP","NTP","TP","SEA","NEA","EA"/)  ;S mean south,N mean North,EA mean East Asia
lats = (/25  ,30   ,35   ,25  , 22  , 30  , 20 /) ;-20
latn = (/30  ,35   ,40   ,40  , 30  , 45  , 45 /)  ;60
lonl = (/70  ,70   ,70   ,70  ,107  ,107  ,107 /)  ;20
lonr = (/105 ,105  ,105  ,105 ,130  ,130  ,122 /);130 ;220
lev  = (/850,500,200/)
var_name = (/"U","V","T"/)
case     = (/"CTRL","NG24h","NG6h"/)
years    = ispan(1979,2005,1)   ;the calculate time
na = 3
nv = 2
nl = 2

path    = "/home/ys17-19/renql"
filein  = path +(/"/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/model/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/model/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
fileout = path + "/project/TP_NUDG/wavelet-20180410/mdata/" + \
            (/"CTRL","NUDG6h","NUDG24h"/) + "_clim_daily_wavelet_" + area(na) + lev(nl) + var_name(nv) + ".nc"

f   = addfile(filein(0)+years(0)+".daily."+var_name(nv)+".nc","r")
time        = cd_calendar(f->time,0)
select_time = ind(time(:,1).ge.1.and.time(:,1).le.12)
ntime       = dimsizes(select_time)
ncase = dimsizes(case)
nyear = dimsizes(years)
var  = new((/nyear,ntime/),float)

;--the parameter for wavelet----------------
mother  = 0    ;0 = Morlet 1 = Paul 2 = DOG (derivative of Gaussian)
param   = 6.0  ;The mother wavelet parameter. If param < 0, then the default is used,
;Morlet k0 (wavenumber), default is 6. Paul m (order), default is 4. DOG m (m-th derivative), default is 2.

N       = ntime
dt      = 1    ;timesteps of time series, the unit here is day
s0      = dt   ;The smallest scale of the wavelet, which is typically equal to 2*dt. s0 = dt for Morlet; s0 = dt/4 for Paul
dj      = 0.25 ;The spacing between discrete scales, which is typically equal to 0.25. A smaller value will give better scale resolution, but will be slower.
jtot    = 1+floattointeger(((log10(N*dt/s0))/dj)/log10(2.))  
;The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj].

npad    = N    ;The total number of points (including padding) to use for the wavelet transform
nadof   = 0    ;Currently ignored (set to zero).

noise   = 1    ;0 means use a white noise for significance test. 1 means use a red noise for significance test. (Most commonly, noise = 1.)
siglvl  = .05  ;The significance level to use
isigtest= 0    ; 0 means do a regular chi-square test, i.e. Eqn (18) from Torrence and Compo. 1 means do a "time-average" test on the global wavelet spectrum.

power  = new((/jtot*npad/),float)
phase  = new((/jtot*npad/),float)
signif = new((/jtot/),float)
gws    = new((/jtot/),float)
mean   = new((/1/),float)
stdev  = new((/1/),float)

do nc = 1, ncase-1,1
print("start calc " + case(nc) + var_name(nv) + lev(nl))
;--------------------------------------------------------------
;read the data
;--------------------------------------------------------------
do nt = 0, nyear-1 ,1
    f   = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    var1 = f->$var_name(nv)$(select_time,{lev(nl)},{lats(na):latn(na)},{lonl(na):lonr(na)})
    var(nt,:) = wgt_areaave_Wrap(var1,1.0,1.0,0) ;get var(nvar,ncase,nyear,365,nlev)
end do
var = var - conform(var,dim_avg_n(var,0),1)  ;remove annual cycle
;var = dim_standardize_n(var,1)

power  = 0
phase  = 0
signif = 0
gws    = 0
mean   = 0
stdev  = 0
do nt = 0, nyear-1 ,1
    wave = wavelet(var(nt,:),mother,dt,param,s0,dj,jtot,npad,noise,isigtest,siglvl,nadof)
    power  = power  + wave@power 
    phase  = phase  + wave@phase 
    signif = signif + wave@signif
    gws    = gws    + wave@gws  
    mean   = mean   + wave@mean
    stdev  = stdev  + wave@stdev
end do
power  = power /nyear
phase  = phase /nyear
signif = signif/nyear
gws    = gws   /nyear
mean   = mean  /nyear
stdev  = stdev /nyear
wave@power  = power 
wave@phase  = phase 
wave@signif = signif
wave@gws    = gws  
wave@mean   = mean 
wave@stdev  = stdev
printVarSummary(wave)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(nc))
ncdf = addfile(fileout(nc),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "the climate wavelet analysis of daily rainfall in " + area(na)
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"part" ,"jtot","time"/)
dimSizes = (/  2    ,jtot  ,N     /)
dimUnlim = (/False  ,False ,False /)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"part","integer"   ,(/"part"/))
filevardef(ncdf,"jtot","integer"   ,(/"jtot"/))
filevardef(ncdf,"time","integer"   ,(/"time"/))
filevardef(ncdf,"wave",typeof(wave) ,(/"part","jtot","time"/))

;Define var Attribute
part    = (/0,1/)
part@long_nmae   = "the real (0,:,:) and imaginary parts (1,:,:) of the wavelet transform"
jtot@long_name   = "The integer number of scales. Scales range from s0 up to s0*2^[(jtot-1)*dj]."
filevarattdef(ncdf,"part",part)
filevarattdef(ncdf,"jtot",jtot)
;filevarattdef(ncdf,"time",time)
filevarattdef(ncdf,"wave",wave)

ncdf->part   = (/0,1/)  
ncdf->jtot   = ispan(1,jtot,1) 
ncdf->time   = ispan(1,N,1)
ncdf->wave   = (/wave/)
end do
end

