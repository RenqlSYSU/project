;******************************************************************************
;use 6 kinds of data to calculate the eof of three jets in winter
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;use the function of eofunc_n_Wrap, eofunc_ts_n_Wrap to 
;calculate the EOF pattern and time serises.
;                                            by Ql_Ren
;                                           2020/02/23
;******************************************************************************
begin
var_name = (/"uwnd","vwnd","omega","air","hgt","ps"/)
nv = 0
opt = True
opt@jopt = 0 ;use covariance matrix (0) or correlation matrix (1) to compute EOFs
neval = 4

area = (/"ME","EA","NA"/);,"GB"
lats = (/ 15 , 20 , 25 /);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/ 40 , 45 , 50 /);, 60  ;90 , 90 , 90 , 90 
lonl = (/ 20 ,105 , 250/);, 0   ;30 ;20
lonr = (/ 60 ,175 , 320/);, 360 ;150;220
narea = dimsizes(area)

path     = "/home/ys17-19/renql/project/2020MiddleEastJet"
filein   = path + "/data/ERA-Interim_DJF_uwnd-year.nc"
fileout  = path + "/data/eof_u_6kinds-" + area + ".nc"

f2 = addfile(path + "/data/ERA-Interim_DJF_ps-year.nc","r")
f  = addfile(filein,"r")
years  = ispan(1979,2015,1)
nyear  = dimsizes(f->season_ave(:,0,0,0))
eof_ts  = new((/6,neval,nyear/),float)
ts_mean = new((/6,neval/),float)

do na = 0,narea-1,1 
vars:= f->season_ave(:,:,{lats(na):latn(na)},{lonl(na):lonr(na)})
printVarSummary(vars)
nlev = dimsizes(vars&lev)
nlat = dimsizes(vars&lat)
nlon = dimsizes(vars&lon)
remove_mean = True
return_info = True
var := dtrend_msg_n(years,vars,remove_mean,return_info,0)
print(area(na) + " : " + lonl(na) + "-" + lonr(na) + "E " + lats(na) + "-" + latn(na) + "N ")
print("max slope is " + max(var@slope))
print("min slope is " + min(var@slope))
print("max intercept is " + max(var@y_intercept))
print("min intercept is " + min(var@y_intercept))
copy_VarMeta(vars,var)
pi   = atan(1.0)*4
wgt  = sqrt(cos( vars&lat*pi/180.0 ))
var  = var*conform(var,wgt,2)
;---------------------------------------------------------------
;different kinds of EOF
;--------------------------------------------------------------
nk = 0
print("")
print("0.horizontal distribution of 200hPa level")
eof0     := eofunc_n(var(:,{200},:,:),neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof0@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var(:,{200},:,:),eof0,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof0@eval(0) + ", " + eof0@eval(1) \
        + ", " + eof0@eval(2) + ", " + eof0@eval(3)) 

nk = 1
print("")
print("1.horizontal distribution of vertical average") 
ptop = 70 
iopt = 0 ;Set to zero. Currently not used.
ps  := f2->season_ave(:,{lats(na):latn(na)},{lonl(na):lonr(na)})
dp  := dpres_plevel(var&lev*100,ps,ptop,iopt)
var1:= dim_sum_n(var*dp,1)/dim_sum_n(dp,1);vertical average
eof1     := eofunc_n(var1,neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof1@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var1,eof1,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof1@eval(0) + ", " + eof1@eval(1) \ 
        + ", " + eof1@eval(2) + ", " + eof1@eval(3)) 

nk = 2
print("")
print("2.vertical section of zonal mean ")
var2 := dim_avg_n(var,3)
eof2     := eofunc_n(var2,neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof2@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var2,eof2,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof2@eval(0) + ", " + eof2@eval(1) \ 
        + ", " + eof2@eval(2) + ", " + eof2@eval(3)) 

nk = 3
print("")
print("3.latitude distribution of vertical average and zonal mean")
var3 := dim_avg_n(var1,2)
eof3     := eofunc_n(var3,neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof3@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var3,eof3,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof3@eval(0) + ", " + eof3@eval(1) \ 
        + ", " + eof3@eval(2) + ", " + eof3@eval(3)) 

nk = 4
print("")
print("4.vertical section of meridional mean") 
var4 := dim_avg_n(var,2)
eof4     := eofunc_n(var4,neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof4@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var4,eof4,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof4@eval(0) + ", " + eof4@eval(1)\ 
        + ", " + eof4@eval(2) + ", " + eof4@eval(3)) 

nk = 5
print("")
print("5.longtitude distribution of vertical average and meridional mean")
var5 := dim_avg_n(var1,1)
eof5     := eofunc_n(var5,neval,opt,0) ;get the (neval,nlat,nlon)
sig      := eofunc_north(eof5@pcvar,nyear,True)
eof_ts_temp   = eofunc_ts_n(var5,eof5,False,0) ;get the (neval,nyear)
eof_ts(nk,:,:)= eof_ts_temp
ts_mean(nk,:) = eof_ts_temp@ts_mean
print("eval : " + eof5@eval(0) + ", " + eof5@eval(1)\ 
        + ", " + eof5@eval(2) + ", " + eof5@eval(3)) 

print("")
do ne = 0,neval-1,1
ne2 = ne + 1
print(area(na) + " EOF" + ne2)
do nk1 = 0    , 4, 1
do nk2 = nk1+1, 5, 1
    r = escorc(eof_ts(nk1,ne,:),eof_ts(nk2,ne,:))
    ;print("correlation between EOF" + nk1 + " and EOF" + nk2 + " neval = " + ne + " : " + r)
    print(" " + r)
end do
print("")
end do
end do

;========================================================================
;save the data
;=======================================================================
system("rm -f " + fileout(na))
ncdf = addfile(fileout(na),"c")
setfileoption(ncdf,"DefineMode",True)

;Set All field attribution
fileAtt =  True
fileAtt@creation_date = systemfunc("date")
fileAtt@discription   = "there are six types of EOF." 
fileattdef(ncdf,fileAtt)  

;Define all Coordinate
dimNames = (/"nkind","neval","year","lev","lat","lon"/)
dimSizes = (/  6    ,neval  ,nyear ,nlev ,nlat ,nlon /)
dimUnlim = (/False  ,False  ,False ,False,False,False/)
filedimdef(ncdf, dimNames, dimSizes, dimUnlim)

;Define var_names, type(string) and dims_names(string)
filevardef(ncdf,"neval","integer",(/"neval"/))
filevardef(ncdf,"year","integer",(/"year"/))
filevardef(ncdf,"lev",typeof(vars&lev),(/"lev"/))
filevardef(ncdf,"lat",typeof(vars&lat),getvardims(vars&lat))
filevardef(ncdf,"lon",typeof(vars&lon),getvardims(vars&lon))
filevardef(ncdf,"eof0",typeof(eof0),(/"neval","lat","lon"/))
filevardef(ncdf,"eof1",typeof(eof1),(/"neval","lat","lon"/))
filevardef(ncdf,"eof2",typeof(eof2),(/"neval","lev","lat"/))
filevardef(ncdf,"eof3",typeof(eof3),(/"neval","lat"/))
filevardef(ncdf,"eof4",typeof(eof4),(/"neval","lev","lon"/))
filevardef(ncdf,"eof5",typeof(eof5),(/"neval","lon"/))
filevardef(ncdf,"eof_ts",typeof(eof_ts),(/"nkind","neval","year"/))
filevardef(ncdf,"ts_mean",typeof(ts_mean),(/"nkind","neval"/))

;Define var Attribute
filevarattdef(ncdf,"lev",vars&lev)
filevarattdef(ncdf,"lat",vars&lat)
filevarattdef(ncdf,"lon",vars&lon)
filevarattdef(ncdf,"eof0",eof0)
filevarattdef(ncdf,"eof1",eof1)
filevarattdef(ncdf,"eof2",eof2)
filevarattdef(ncdf,"eof3",eof3)
filevarattdef(ncdf,"eof4",eof4)
filevarattdef(ncdf,"eof5",eof5)
filevarattdef(ncdf,"eof_ts",eof_ts)
filevarattdef(ncdf,"ts_mean",ts_mean)

ncdf->lev    = (/vars&lev/) 
ncdf->lat    = (/vars&lat/)
ncdf->lon    = (/vars&lon/)
ncdf->eof0 = (/eof0/)
ncdf->eof1 = (/eof1/)
ncdf->eof2 = (/eof2/)
ncdf->eof3 = (/eof3/)
ncdf->eof4 = (/eof4/)
ncdf->eof5 = (/eof5/)
ncdf->eof_ts = (/eof_ts/)
ncdf->ts_mean = (/ts_mean/)
end do
end

