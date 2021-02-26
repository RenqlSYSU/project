;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path     = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
filein   = path + "ERA-Interim_DJF_sfvp-year.nc" 
fileout  = path + "ERA-Interim_DJF_rws-reg2.nc"
var_name = (/"ud","vd","vort"/)

g  = 9.8 ;m/(s*s)
cp = 1004.0 ; J/(K kg)  [ m2/(K s2) ] 
R  = 287.0  ; J/(kg*K) for dry air
a  = 6378388 ;the radius of earth, m
sclhgt = 7991 ;m 
pi = atan(1.0)*4
print(pi)
years = ispan(1979,2016,1)
nyear = dimsizes(years)-1
years@calendar = "standard"

;area: east Asia
lats =  5  ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/850,500,200,100/)

;==========================================================
;read data
;============================================================
;------------------EOF time series------------------------
neval = 3
area  = (/"ME","EA","NA"/);,"GB"
narea = dimsizes(area)
ts    = new((/narea,neval,nyear/),float)
eval  = new((/narea,neval/),float)
pcvar = new((/narea,neval/),float)
sig   = new((/narea,neval/),logical)
do na = 0, narea-1,1
    f  = addfile(path + "eof_u_6kinds-" + area(na) + ".nc","r")
    eof0 := f->eof0
    eval(na,:)  = eof0@eval(0:(neval-1))
    pcvar(na,:) = eof0@pcvar(0:(neval-1))
    ts(na,:,:)  = f->eof_ts(0, 0:(neval-1), 0:(nyear-1))
    sig(na,:)   = eofunc_north(pcvar(na,:),nyear,False)
end do
ts = - ts/conform(ts,sqrt(eval),(/0,1/))
delete(eof0)

;--------------------------------------------------------
f    = addfile(filein,"r")
ud   = f->$var_name(0)$(:,{lev},{lats:latn},{lonl:lonr})
vd   = f->$var_name(1)$(:,{lev},{lats:latn},{lonl:lonr})
vort = f->$var_name(2)$(:,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(ud&lat)
nlon  = dimsizes(vd&lon)
nlev  = dimsizes(lev)
print(nyear)

iyear = 0
ilev  = 1
ilat  = 2
ilon  = 3
lat  = ud&lat
lon  = ud&lon
lat  = lat*pi/180.0
lon  = lon*pi/180.0
opt    = 0     ;used by center_finite_diff_n, no meanging 
loncyc = True ;used by center_finite_diff_n
cyclic = False ;used by center_finite_diff_n
;means NOT treated as cyclic in differenced coordinates
;the end values will use a one-sided difference scheme for the end points. 

ano_ud   = regCoef_n(ts, ud,  2, iyear) ;narea,neval,nlev,nlat,nlon
ano_vd   = regCoef_n(ts, vd,  2, iyear) ;narea,neval,nlev,nlat,nlon
clm_vort = conform(ano_ud,dim_avg_n(vort,iyear),(/ilev+1,ilat+1,ilon+1/))
f0 = conform(ano_ud, 2.0*(2.0*pi/24.0/3600.0)*sin(lat), ilat+1)
coslat = conform(ano_ud, cos(lat), ilat+1)

term1 = (f0+clm_vort)*(center_finite_diff_n(ano_ud,lon,loncyc,opt,ilon+1) \
        + center_finite_diff_n(ano_vd*coslat,lat,cyclic,opt,ilat+1))/(a*coslat) 
term1 = (ano_ud*center_finite_diff_n(f0+clm_vort,lon,loncyc,opt,ilon+1) \
        + ano_vd*center_finite_diff_n((f0+clm_vort)*coslat,lat,cyclic,opt,ilat+1))/(a*coslat) + term1
delete([/f0, clm_vort, ano_vd, ano_ud/])

ano_vort = regCoef_n(ts, vort,2, iyear) ;narea,neval,nlev,nlat,nlon
clm_ud   = conform(ano_vort,dim_avg_n(ud,iyear),(/ilev+1,ilat+1,ilon+1/))
clm_vd   = conform(ano_vort,dim_avg_n(vd,iyear),(/ilev+1,ilat+1,ilon+1/))
delete([/vd,vort/])

term1 = ano_vort*(center_finite_diff_n(clm_ud,lon,loncyc,opt,ilon+1) \
        + center_finite_diff_n(clm_vd*coslat,lat,cyclic,opt,ilat+1))/(a*coslat)  + term1
term1 = (clm_ud*center_finite_diff_n(ano_vort,lon,loncyc,opt,ilon+1) \
        + clm_vd*center_finite_diff_n(ano_vort*coslat,lat,cyclic,opt,ilat+1))/(a*coslat) + term1
delete([/ano_vort, clm_vd, clm_ud, coslat/])
rws   = -term1
delete(term1)

rws@long_name = "Rossby wave source"
rws@units   = "s^-2" 
rws@pcvar   = ndtooned(pcvar)
rws@eval    = ndtooned(eval)
rws@sig     = ndtooned(sig)

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout)
setfileoption("nc","Format","LargeFile") 
ncdf = addfile(fileout,"c")
setfileoption(ncdf,"definemode",True)

;set all field attribute
fileatt =  True
fileatt@creation_date = systemfunc("date")
fileatt@discription   = "Rossby wave source, ts*-1"
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"area","neval","lev","lat","lon"/);
dimSizes = (/narea ,neval  ,nlev ,nlat ,nlon/) ;
dimUnlim = (/False ,False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(ud&lev),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(ud&lat),getvardims(ud&lat))
filevardef(ncdf,"lon" ,typeof(ud&lon),getvardims(ud&lon))
filevardef(ncdf,"rws",typeof(rws),(/"area","neval","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",ud&lev)
filevarattdef(ncdf,"lat",ud&lat)
filevarattdef(ncdf,"lon",ud&lon)
filevarattdef(ncdf,"rws",rws)

ncdf->lev  = (/lev/)
ncdf->lat  = (/ud&lat/)
ncdf->lon  = (/ud&lon/)
ncdf->rws  = (/rws/)
end

