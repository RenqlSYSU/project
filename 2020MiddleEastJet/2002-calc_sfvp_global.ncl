;***************************************************************************
;use the reanalysis data to calculate diabatic heating, transient eddy heating
;and Transient eddy vorticity forcing 
;                                            by Ql_Ren
;                                           2020/02/04
;******************************************************************************
begin
path1   = "/home/ys17-19/data-observation/ERA-interim/pressure/"
path2   = "/home/ys17-19/renql/project/2020MiddleEastJet/data/"
path3   = "/home/ys17-19/renql/model/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.1979.daily.U.nc" ;to get time attribute
var_name = (/"uwnd","vwnd"/)

years    = ispan(1979,2016,1)
years@calendar = "standard"

;area: east Asia
lats = -90 ;-20
latn =  90 ;60
lonl =  0  ;20
lonr = 360 ;220
lev  = (/1000,925,850,700,600, 500,400,350,300,250, 200,175,150,125,100/)

;==========================================================
;create an array(nvar,ncase,nyear) to store the data
;============================================================
f_model = addfile(path3,"r")
time = cd_calendar(f_model->time,0)
delete(f_model)
printVarSummary(time)

f    = addfile(path1+var_name(0)+"/"+var_name(0)+".interim."+years(0)+".nc", "r")
vars = f->$var_name(0)$(0,{lev},{lats:latn},{lonl:lonr})
nlat  = dimsizes(vars&lat)
nlon  = dimsizes(vars&lon)
nlev  = dimsizes(lev)
nyear = dimsizes(years)

select_time  = ind(time(:,1).eq.12.or.time(:,1).le.2)
ntime = dimsizes(select_time)
select_time1 = select_time
select_time1(58:89) = select_time(58:89) + 1

year_option = True;False ;
fileout  = path2 + "ERA-Interim_DJF_sfvp_global-year.nc"
;==========================================================
;read the data
;============================================================
u   = new((/nyear,ntime,nlev,nlat,nlon/),float)
v   = new((/nyear,ntime,nlev,nlat,nlon/),float)
do nt = 0,nyear-1,1
print("read: " + years(nt))
if(isleapyear(years(nt))) then ;leap year,366 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time1,{lev},{lats:latn},{lonl:lonr})
else   ;non-leap year, 365 day
    nv = 0
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    u(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
    nv = 1
    f      = addfile(path1+var_name(nv)+"/"+var_name(nv)+".interim."+years(nt)+".nc","r")
    v(nt,:,:,:,:)  = f->$var_name(nv)$(select_time,{lev},{lats:latn},{lonl:lonr})
end if
end do

;--------------------------------------
;DJF
;---------------------------------------------------
u2   := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
v2   := new((/nyear-1,ntime,nlev,nlat,nlon/),float)
do ny = 0,nyear-2,1
    u2(ny,0:30,:,:,:)  = u(ny ,59:89,:,:,:)
    u2(ny,31:89,:,:,:) = u(ny+1,0:58,:,:,:)
    v2(ny,0:30,:,:,:)  = v(ny ,59:89,:,:,:)
    v2(ny,31:89,:,:,:) = v(ny+1,0:58,:,:,:)
end do
u3 = dim_avg_n(u2,1)
v3 = dim_avg_n(v2,1)
delete([/u,v,u2,v2/])
printVarSummary(v3) ;(/nyear-1,nlev,nlat,nlon/)

if(any(ismissing(u3))) then
    print("uwnd has missing")
end if
if(any(ismissing(v3))) then
    print("vwnd has missing")
end if

vort = new((/nyear-1,nlev,nlat,nlon/),float)
divg = new((/nyear-1,nlev,nlat,nlon/),float)
sf   = new((/nyear-1,nlev,nlat,nlon/),float)
vp   = new((/nyear-1,nlev,nlat,nlon/),float)
ur   = new((/nyear-1,nlev,nlat,nlon/),float)
vr   = new((/nyear-1,nlev,nlat,nlon/),float)
ud   = new((/nyear-1,nlev,nlat,nlon/),float)
vd   = new((/nyear-1,nlev,nlat,nlon/),float)
uv2vrdvf(u3,v3,vort,divg) ;Computes the vorticity and divergence via spherical harmonics
uv2sfvpf(u3,v3,sf,vp) ;calculate stream function and velocity potential via spherical harmonics
vr2uvf(vort,ur,vr)    ;calculate the rotational wind components via spherical harmonics
dv2uvf(divg,ud,vd)    ;calculate the divergent wind components via spherical harmonics

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
fileatt@discription   = "DJF mean, all variables are calculated by spherical harmonics" 
fileattdef(ncdf,fileatt)  

;define all coordinate
dimNames = (/"year" ,"lev","lat","lon"/);
dimSizes = (/nyear-1,nlev ,nlat ,nlon/) ;
dimUnlim = (/False  ,False,False,False/);
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;define var_names, type(string) and dims_names(string)
filevardef(ncdf,"lev" ,typeof(vars&level),(/"lev"/))
filevardef(ncdf,"lat" ,typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon" ,typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"vort",typeof(vort),(/"year","lev","lat","lon"/))
filevardef(ncdf,"divg",typeof(divg),(/"year","lev","lat","lon"/))
filevardef(ncdf,"sf",typeof(sf),(/"year","lev","lat","lon"/))
filevardef(ncdf,"vp",typeof(vp),(/"year","lev","lat","lon"/))
filevardef(ncdf,"ur",typeof(ur),(/"year","lev","lat","lon"/))
filevardef(ncdf,"vr",typeof(vr),(/"year","lev","lat","lon"/))
filevardef(ncdf,"ud",typeof(ud),(/"year","lev","lat","lon"/))
filevardef(ncdf,"vd",typeof(vd),(/"year","lev","lat","lon"/))

;define var attribute
filevarattdef(ncdf,"lev",vars&level)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"vort",vort)
filevarattdef(ncdf,"divg",divg)
filevarattdef(ncdf,"sf",sf)
filevarattdef(ncdf,"vp",vp)
filevarattdef(ncdf,"ur",ur)
filevarattdef(ncdf,"vr",vr)
filevarattdef(ncdf,"ud",ud)
filevarattdef(ncdf,"vd",vd)

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
end

