module init_topo_mod
contains
subroutine init_topo
!**** *init_topo*  - init settings 

!     --------------
!        original : 2011-05-04
!                   based on the routines by dai yamazaki
!     ------------------------------------------------------------------

use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam,   tmpnam,   nx,ny,nlfp, imis
use mod_input  ,only: cgrarea,celevtn,cnxtdst,crivwth,crivlen,crivhgt,cfldhgt,&
                      crivclinc,crivparnc,lfld,lkine,lmapcdf,lmapend

use mod_map    ,only: d2nxtdst, d2grarea, d2elevtn, d2rivlen, d2rivwth, d2rivhgt, d2fldhgt, &
                      d2rivelv, d2fldgrd, d2rivstomax, d2fldstomax, dfrcinc
use mod_map    ,only: i2vector, nseqall,  nseqmax
implicit none

integer(kind=jpim) :: iseq, ix, iy

!!=========================
write(lognam,*) " "
write(lognam,*) "*** init_topo: ***"

allocate( d2grarea(nseqmax,1) )
allocate( d2elevtn(nseqmax,1) )
allocate( d2nxtdst(nseqmax,1) )
allocate( d2rivlen(nseqmax,1) )
allocate( d2rivwth(nseqmax,1) )
allocate( d2rivhgt(nseqmax,1) )
allocate( d2fldhgt(nseqmax,1,nlfp) )

d2grarea(:,:)=0.d0
d2elevtn(:,:)=0.d0
d2nxtdst(:,:)=0.d0
d2rivlen(:,:)=0.d0
d2rivwth(:,:)=0.d0
d2fldhgt(:,:,:)=0.d0

if ( .not. lmapcdf ) then
  call read_topo_bin
else
  call read_topo_cdf
endif

write(lognam,*) 'init_topo: calc rivstomax (flood initiation)'

allocate(d2rivstomax(nseqmax,1))
if ( lfld ) then
  d2rivstomax(:,:) = d2rivlen(:,:) * d2rivwth(:,:) * d2rivhgt(:,:)
else
  d2rivstomax(:,:) = 1.d18
  write(lognam,*) 'init_topo: no floodplain (rivstomax=1.d18)'
endif

write(lognam,*) 'init_topo: calc river bed elevation (rivelv)'
allocate(d2rivelv(nseqmax,1))
d2rivelv(:,:) = d2elevtn(:,:) - d2rivhgt(:,:)

write(lognam,*) 'init_topo: calc floodplain parameters'
allocate(d2fldstomax(nseqmax,1,nlfp))
allocate(d2fldgrd(nseqmax,1,nlfp))
call set_fldstg

contains




!!==================================================
subroutine set_fldstg
implicit none

integer(kind=jpim)  ::  iseq, i
real(kind=jprb)     ::  dstonow
real(kind=jprb)     ::  dstopre
real(kind=jprb)     ::  dhgtpre
real(kind=jprb)     ::  dwthinc
! ================================================
d2fldstomax(:,:,:) = 0.d0
d2fldgrd(:,:,:)    = 0.d0
dfrcinc=dble(nlfp)**(-1.)
!
do iseq=1, nseqall
  dstopre = d2rivstomax(iseq,1)
  dhgtpre = 0.d0
  dwthinc = d2grarea(iseq,1) * d2rivlen(iseq,1)**(-1.) * dfrcinc
  do i=1, nlfp
    dstonow = d2rivlen(iseq,1) * ( d2rivwth(iseq,1) + dwthinc*(dble(i)-0.5) ) * (d2fldhgt(iseq,1,i)-dhgtpre)
    d2fldstomax(iseq,1,i) = dstopre + dstonow
    d2fldgrd(iseq,1,i) = (d2fldhgt(iseq,1,i)-dhgtpre) * dwthinc**(-1.)
    dstopre = d2fldstomax(iseq,1,i)
    dhgtpre = d2fldhgt(iseq,1,i)
  end do
end do
!
return
end subroutine set_fldstg





! =================================================
subroutine read_topo_bin
implicit none

real(kind=jprm),allocatable  :: r2temp(:,:)
real(kind=jprb),allocatable  :: d2temp(:,:)
integer(kind=jpim)           :: ilfp

allocate(r2temp(nx,ny))
allocate(d2temp(nseqmax,1))
!!===================================
!! 2. read binary data

write(lognam,*)'init_topo: unit-catchment area : ',trim(cgrarea) 
open(tmpnam,file=cgrarea,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2grarea)
close(tmpnam)

write(lognam,*)'init_topo: ground elevation : ',trim(celevtn)
open(tmpnam,file=celevtn,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2elevtn)
close(tmpnam)

write(lognam,*)'init_topo: downstream distance : ',trim(cnxtdst)
open(tmpnam,file=cnxtdst,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2nxtdst)
close(tmpnam)

write(lognam,*)'init_topo: river channel width : ',trim(crivwth)
open(tmpnam,file=crivwth,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2rivwth)
close(tmpnam)

write(lognam,*)'init_topo: river channel depth : ',trim(crivhgt)
open(tmpnam,file=crivhgt,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2rivhgt)
close(tmpnam)

write(lognam,*)'init_topo: river channel length : ',trim(crivlen)
open(tmpnam,file=crivlen,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp(:,:)
  if( lmapend ) call conv_end(r2temp)
call map2vec(r2temp,d2rivlen)
close(tmpnam)

write(lognam,*)'init_topo: floodplain elevation profile : ',trim(cfldhgt)
open(tmpnam,file=trim(cfldhgt),form='unformatted',access='direct',recl=4*nx*ny)
do ilfp=1,nlfp
  read(tmpnam,rec=ilfp) r2temp
  if( lmapend ) call conv_end(r2temp)
  call map2vec(r2temp,d2temp)
  d2fldhgt(:,:,ilfp)= d2temp(:,:)
enddo
close(tmpnam)

deallocate(r2temp)
deallocate(d2temp)

end subroutine read_topo_bin



! =================================================
subroutine read_topo_cdf
#ifdef UseCDF
use netcdf 
use lib_netcdf_util
implicit none

integer(kind=jpim)           :: ncid,varid
integer(kind=jpim)           :: ilev

real(kind=jprm),allocatable     :: r2temp(:,:)
real(kind=jprb),allocatable     :: d2temp(:,:)

!!=============================
allocate(r2temp(nx,ny))
allocate(d2temp(nseqmax,1))

!! clim file
call ncerror (nf90_open(crivclinc,nf90_nowrite,ncid),'opening '//trim(crivclinc) )

write(lognam,*)'init_topo: grarea:',trim(crivclinc)
call ncerror ( nf90_inq_varid(ncid,'grarea',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' )
call map2vec(r2temp,d2grarea)

write(lognam,*)'init_topo: elevtn:',trim(crivclinc)
call ncerror ( nf90_inq_varid(ncid,'elevtn',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' ) 
call map2vec(r2temp,d2elevtn)

write(lognam,*)'init_topo: nxtdst:',trim(crivclinc)
call ncerror ( nf90_inq_varid(ncid,'nxtdst',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' ) 
call map2vec(r2temp,d2nxtdst)

write(lognam,*)'init_topo: rivlen:',trim(crivclinc)
call ncerror ( nf90_inq_varid(ncid,'rivlen',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' ) 
call map2vec(r2temp,d2rivlen)

write(lognam,*)'init_topo: fldhgt:',trim(crivclinc)
call ncerror ( nf90_inq_varid(ncid,'fldhgt',varid),'getting id' )
do ilev=1,nlfp
  call ncerror ( nf90_get_var(ncid,varid,r2temp,(/1,1,ilev/),(/nx,ny,1/)),'reading data' ) 
  call map2vec(r2temp,d2temp)
  d2fldhgt(:,:,ilev)=d2temp(:,:)
enddo

call ncerror( nf90_close(ncid))

!!========== 
!! par file 
call ncerror (nf90_open(crivparnc,nf90_nowrite,ncid),'opening '//trim(crivparnc) )

write(lognam,*)'init_topo: rivwth:',trim(crivparnc)
call ncerror ( nf90_inq_varid(ncid,'rivwth',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' ) 
call map2vec(r2temp,d2rivwth)

write(lognam,*)'init_topo: rivhgt:',trim(crivparnc)
call ncerror ( nf90_inq_varid(ncid,'rivhgt',varid),'getting id' )
call ncerror ( nf90_get_var(ncid,varid,r2temp),'reading data' ) 
call map2vec(r2temp,d2rivhgt)

call ncerror( nf90_close(ncid))

deallocate(r2temp)
deallocate(d2temp)

#endif

end subroutine read_topo_cdf




! =================================================
subroutine conv_end(r2temp)

implicit none
real(kind=jprm)    :: r2temp(nx,ny)

#ifdef convend
do iy=1, ny
  do ix=1, ny
    call endian4(r2temp(ix,iy))
  end do
end do
#endif

end subroutine conv_end



! =================================================
subroutine map2vec(r2temp,d2var)
implicit none

real(kind=jprm)              :: r2temp(nx,ny)
real(kind=jprb)              :: d2var(nseqmax,1)

do iy=1, ny
  do ix=1, nx
    if( i2vector(ix,iy)>0 )then
      iseq=i2vector(ix,iy)
      d2var(iseq,1)=dble(r2temp(ix,iy))
    endif
  end do
end do

end subroutine map2vec




end subroutine init_topo
end module init_topo_mod
