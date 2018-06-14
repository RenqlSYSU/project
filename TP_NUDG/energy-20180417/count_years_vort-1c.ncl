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
nc = 2

;-----area of east asia 
lats = 20
latn = 45
lonl = 65
lonr = 125
lev  = (/850,500,200/)
nl = 1
fileout = path + "/energy-20180417/mdata/JJAvort_number" + case(nc) + lev(nl) + ".txt"

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

do nc = 0,ncase-1,1
print(" ")
print("=================================================================================")
print(" the number of vortices generating in " + case(nc) + " " + lev(nl) + "hPa")
print(" ")
NUM_a = 0
NUM_T = 0
NUM_av1 = 0
NUM_TV1 = 0
NUM_av2 = 0
NUM_TV2 = 0
do nt = 0,nyear-1,1
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
    
    do nx = 1,nlon-2,1
    do ny = 1,nlat-2,1
    grads(:,ny,nx) = (hgt(:,ny,nx-1)-hgt(:,ny,nx));/(cos(lat(ny)*3.14/180)*6371393*(lon(nx-1)-lon(nx))) ;<0
    grade(:,ny,nx) = (hgt(:,ny,nx+1)-hgt(:,ny,nx));/(cos(lat(ny)*3.14/180)*6371393*(lon(nx+1)-lon(nx))) ;>0
    gradn(:,ny,nx) = (hgt(:,ny+1,nx)-hgt(:,ny,nx));/(6371393*(lat(ny+1)-lat(ny))) ;>0
    gradw(:,ny,nx) = (hgt(:,ny-1,nx)-hgt(:,ny,nx));/(6371393*(lat(ny-1)-lat(ny))) ;<0
    end do
    end do
    copy_VarMeta(hgt,grads)
    copy_VarMeta(hgt,grade)
    copy_VarMeta(hgt,gradn)
    copy_VarMeta(hgt,gradw)
    
    grads1d = ndtooned(grads)
    grade1d = ndtooned(grade)
    gradn1d = ndtooned(gradn)
    gradw1d = ndtooned(gradw)
    indices := ind_resolve(ind(grads1d.gt.0.and.grade1d.gt.0.and.gradw1d.gt.0.and.gradn1d.gt.0),dimsizes(grads)) ;two dimensions
    dim_id  := dimsizes(indices)
    number  := dim_id(0)
    numberTP:= dimsizes(ind(lat(indices(:,1)).gt.25.and.lat(indices(:,1)).lt.40.and.lon(indices(:,2)).gt.70.and.lon(indices(:,2)).lt.105))
    
    num_vort1   = 0
    num_vortTP1 = 0
    num_vort2   = 0
    num_vortTP2 = 0
    do nn = 0,number-1,1
        itime= indices(nn,0)
        ilat = indices(nn,1)
        ilon = indices(nn,2)
        if(vort(itime,ilat,ilon).gt.0) then
        num_vort1 = num_vort1 + 1
        if(lat(ilat).gt.25.and.lat(ilat).lt.40.and.lon(ilon).gt.70.and.lon(ilon).lt.105) then
        num_vortTP1 = num_vortTP1 + 1
        end if
        end if
        
        if(vort(itime,ilat,ilon).gt.0.and.vort(itime,ilat+1,ilon).gt.0.and.vort(itime,ilat-1,ilon).gt.0 \
            .and.vort(itime,ilat,ilon-1).gt.0.and.vort(itime,ilat,ilon+1).gt.0) then
        num_vort2 = num_vort2 + 1
        if(lat(ilat).gt.25.and.lat(ilat).lt.40.and.lon(ilon).gt.70.and.lon(ilon).lt.105) then
        num_vortTP2 = num_vortTP2 + 1
        end if
        end if
    end do
    print(years(nt) + " " + number + " " + numberTP + " " + num_vort1 + " " + num_vortTP1 + " " + num_vort2 + " " + num_vortTP2)
    
    NUM_a = NUM_a + number
    NUM_T = NUM_T + numberTP
    NUM_av1 = NUM_av1 + num_vort1
    NUM_av2 = NUM_av2 + num_vort2
    NUM_TV1 = NUM_TV1 + num_vortTP1
    NUM_TV2 = NUM_TV2 + num_vortTP2
end do
print(" ")
print("all year in east Asia: " + NUM_a + ",all year in TP: " + NUM_T + " " + NUM_av1 + " " + NUM_TV1 + " " + NUM_av2 + " " + NUM_TV2)
print("=================================================================================")
print(" ")
end do
end
