;***************************************************************************
;Draw the first three EOF  principal components series 
;   1.horizontal distribution of 200hPa level
;   2.horizontal distribution of vertical average 
;   3.vertical section of zonal mean 
;   4.latitude distribution of vertical average and zonal mean
;   5.vertical section of meridional mean 
;   6.longtitude distribution of vertical average and meridional mean
;                                            by Ql_Ren
;                                           2020/06/12
;******************************************************************************
begin
area = (/"ME","EA" ,"NA" /);,"GB"
lats = (/15.0, 20.0, 25.0/);, 10  ;-90 ,-90 ,-90 ,-90 
latn = (/40.0, 45.0, 50.0/);, 60  ;90 , 90 , 90 , 90 
lonl = (/20.0,105.0,250.0/);, 0   ;30 ;20
lonr = (/60.0,175.0,320.0/);, 360 ;150;220
years = ispan(1979,2015,1)
nyear = dimsizes(years)
narea = dimsizes(area)
neval = 4
nkind = 6
path   = "/home/ys17-19/renql/project/2020MiddleEastJet"
filename   = path + "/data/eof_u_6kinds-" + area + ".nc"

eof_ts = new((/narea, nkind, neval, nyear/),float)

do na = 0, narea-1, 1
    f = addfile(filename(na),"r")
    eof_ts(na,:,:,:) = f->eof_ts
    lat := f->lat
    lon := f->lon
    nlat = dimsizes(lat)
    nlon = dimsizes(lon)
    lats(na) = lat(0)
    latn(na) = lat(nlat-1)
    lonl(na) = lon(0)
    lonr(na) = lon(nlon-1)
print(area(na) + " : " + lonl(na) + "-" + lonr(na) + "E " + lats(na) + "-" + latn(na) + "N ")
end do
print(" ")
print(" ")

do nk = 0,nkind-1,1
do na = 0,1,1;narea-1, 1
    do ne = 0,2,1 ;
    do na2 = na+1,narea-1,1
    do ne2 = 0,2,1
    r = escorc(eof_ts(na,nk,ne,:),eof_ts(na2,nk,ne2,:))
    print("type" + nk + " correlation between " + area(na) + " EOF" + ne + " and " + area(na2) + " EOF" + ne2 + " = " + r)
    end do
    end do
    print(" ")
    end do
end do 
print(" ")
print(" ")
end do
end

