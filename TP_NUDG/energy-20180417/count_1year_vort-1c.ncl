;******************************************************************************
;mission:
;1 read the U,V,Z3 and time over East Asia every year
;2 calculate the vorticity and the gradient of pressure
;3 Judging the low pressure area and count the plateau vortices and date
;4 Judging the eastward vortices and count the number and data
;                                            by Ql_Ren
;                                           2017/12/19
;******************************************************************************
begin
;------the path of the file
path    = "/users/yangsong3/renql/project/TP_NUDG"
filein  = path +(/"/data/AMIP-CTRL/AMIP_C5PM.cam.h1.YEAR.",\
                  "/data/TP-NUDG-24h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR.",\
                  "/data/TP-NUDG-6h/AMIP_C5PM_TP_NUDG.cam.h1.YEAR."/)
years    = ispan(1979,2005,1)   ;the calculate time
case     = (/"CTRL","NG24h","NG6h"/)
var_name = (/"U","V","Z3"/)
nt = 19
nc = 2
fileout = path + "/energy-20180417/mdata/JJAvort_number" + case(nc) + years(nt) + ".txt"

;-----area of east asia 
lats = 20
latn = 45
lonl = 65
lonr = 125
lev  = (/850,500,200/)
nl = 1

;-----area of TP 
maxlon = 105
minlon = 70
maxlat = 40
minlat = 25

;----find the time index of the JJA and define the dimension of the var
f          = addfile(filein(0) + years(0) + ".daily." + var_name(0)+".nc","r")
time       = cd_calendar(f->time,0)
select_time= ind(time(:,1).ge.6.and.time(:,1).le.8)
ntime      = dimsizes(select_time)
lat   = f->lat({lats:latn})
lon   = f->lon({lonl:lonr})
nlat  = dimsizes(lat)
nlon  = dimsizes(lon)
nyear = dimsizes(years)
ncase = dimsizes(case)
grads  = new((/ntime,nlat,nlon/),float) ;south
grade  = new((/ntime,nlat,nlon/),float) ;east
gradn  = new((/ntime,nlat,nlon/),float) ;north
gradw  = new((/ntime,nlat,nlon/),float) ;west

;do nc = 0,ncase-1,1
;do nt = 0,nyear-1,1
    nv = 0
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    u = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})
    nv = 1
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    v = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})
    nv = 2
    f = addfile(filein(nc)+years(nt)+".daily."+var_name(nv)+".nc","r")
    hgt = f->$var_name(nv)$(select_time,{lev(nl)},{lats:latn},{lonl:lonr})

    vort = uv2vr_cfd(u, v, v&lat, v&lon, 2)*1000000 ;calculate the vorticity
    copy_VarMeta(u,vort)
    printVarSummary(hgt)
    printVarSummary(u)
    printVarSummary(v)
    printVarSummary(lat)
    
    do nx = 1,nlon-2,1
    do ny = 1,nlat-2,1
    grads(:,ny,nx) = (hgt(:,ny,nx-1)-hgt(:,ny,nx));/(cos(lat(ny)*3.14/180)*6371393*(lon(nx-1)-lon(nx))) ;<0
    grade(:,ny,nx) = (hgt(:,ny,nx+1)-hgt(:,ny,nx));/(cos(lat(ny)*3.14/180)*6371393*(lon(nx+1)-lon(nx))) ;>0
    gradn(:,ny,nx) = (hgt(:,ny+1,nx)-hgt(:,ny,nx));/(6371393*(lat(ny+1)-lat(ny))) ;>0
    gradw(:,ny,nx) = (hgt(:,ny-1,nx)-hgt(:,ny,nx));/(6371393*(lat(ny-1)-lat(ny))) ;<0
    end do
    end do
    copy_VarMeta(u,grads)
    copy_VarMeta(u,grade)
    copy_VarMeta(u,gradn)
    copy_VarMeta(u,gradw)
    
    grads1d = ndtooned(grads)
    grade1d = ndtooned(grade)
    gradn1d = ndtooned(gradn)
    gradw1d = ndtooned(gradw)
    indices = ind_resolve(ind(grads1d.gt.0.and.grade1d.gt.0.and.gradw1d.gt.0.and.gradn1d.gt.0),dimsizes(grads)) ;two dimensions
    dim_id  = dimsizes(indices)
    number  = dim_id(0)
    print(years(nt) + "YEAR " + number + " vortices generate in " + case(nc) + lev(nl) )
    print("===============================================================")

    numberTP= dimsizes(ind(lat(indices(:,1)).gt.25.and.lat(indices(:,1)).lt.40.and.lon(indices(:,2)).gt.70.and.lon(indices(:,2)).lt.105))
    print(years(nt) + "YEAR " + numberTP + "Plateau vortices generate" + case(nc) + lev(nl) )
    print("===============================================================")

    numberTP = 0
    do nn = 0,number-1,1
    if(lat(indices(nn,1)).gt.25.and.lat(indices(nn,1)).lt.40.and.lon(indices(nn,2)).gt.70.and.lon(indices(nn,2)).lt.105) then
    numberTP = numberTP + 1
    print(years(nt) + "/" + time(select_time(indices(nn,0)),1) + "/" + time(select_time(indices(nn,0)),2) \
            + " " + hgt(indices(nn,0),indices(nn,1),indices(nn,2)) + "m " + lon(indices(nn,2)) + "E " + lat(indices(nn,1)) + "N")
    end if
    end do
    print("===============================================================")
    print(years(nt) + "YEAR " + numberTP + " Plateau vortices generate" + case(nc) + lev(nl) )
;end do
;end do
end
