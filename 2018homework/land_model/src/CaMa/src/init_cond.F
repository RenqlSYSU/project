module init_cond_mod
contains 
subroutine init_cond(iyyyy, imm, idd)
! ===============================================
! set initial condition
! ===============================================
use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam,nx,ny,nlfp,dmis,&
                      irestart,crestdir,creststo,lrestcdf,lstoonly,lpthout
use mod_map    ,only: nseqmax, npthout, npthlev, i2vector
use mod_prog   ,only: d2rivsto_in, d2fldsto_in, d2rivout_in, d2fldout_in, d1pthflw_in
use mod_prog   ,only: d2rivsto_out,d2fldsto_out,d2rivout_out,d2fldout_out,d1pthflw_out, d2runoff, &
                      d2rivout_pre,d2fldout_pre,d1pthflw_pre,d2rivdph_pre,d2fldsto_pre
use mod_output ,only: cvnrestart,csufbin,csufcdf
#ifdef UseMPI
use mod_map    ,only: regionthis
#endif

use calc_fldstg_mod

implicit none

integer(kind=jpim), intent(in) :: iyyyy, imm, idd   !! start date

write(lognam,*) " "
write(lognam,*) "****init_cond: ****"

!!=============
!! 1. allocate 
allocate(d2runoff(nseqmax,1))
allocate(d2rivsto_in(nseqmax,1))
allocate(d2fldsto_in(nseqmax,1))

allocate(d2rivout_in(nseqmax,1))
allocate(d2fldout_in(nseqmax,1))


allocate(d2rivsto_out(nseqmax,1))
allocate(d2fldsto_out(nseqmax,1))

allocate(d2rivout_out(nseqmax,1))
allocate(d2fldout_out(nseqmax,1))

allocate(d2rivout_pre(nseqmax,1))
allocate(d2fldout_pre(nseqmax,1))
allocate(d2rivdph_pre(nseqmax,1))
allocate(d2fldsto_pre(nseqmax,1))

allocate(d1pthflw_in(npthout,npthlev))
allocate(d1pthflw_out(npthout,npthlev))
allocate(d1pthflw_pre(npthout,npthlev))
!!========================
!! 2. restart or not ! 
if ( irestart .eq. 1 ) then
  if ( lrestcdf ) then
    call read_rest_cdf
  else
    call read_rest_bin
  endif

  if( lstoonly )then
    d2rivout_out(:,:)=0._jprb
    d2fldout_out(:,:)=0._jprb

    d2rivout_pre(:,:)=0._jprb
    d2fldout_pre(:,:)=0._jprb

    d2fldsto_pre(:,:)=d2fldsto_out(:,:)

    d1pthflw_out(:,:)=0._jprb
    d1pthflw_pre(:,:)=0._jprb

    d2rivdph_pre(:,:)=0._jprb     !! bugfix v362
  endif

else
  ! set to zero ...
  write(lognam,*)'read_rest: no restart... setting to zero'
  d2rivsto_out(:,:)=0._jprb
  d2fldsto_out(:,:)=0._jprb

  d2rivout_out(:,:)=0._jprb
  d2fldout_out(:,:)=0._jprb

  d2rivout_pre(:,:)=0._jprb
  d2fldout_pre(:,:)=0._jprb

  d2rivdph_pre(:,:)=0._jprb
  d2fldsto_pre(:,:)=0._jprb

  d1pthflw_out(:,:)=0._jprb
  d1pthflw_pre(:,:)=0._jprb

  ! set initial water surface elevation to sea surface level
  call storage_sea_surface

endif

d2rivsto_in(:,:)=d2rivsto_out(:,:)
d2fldsto_in(:,:)=d2fldsto_out(:,:)
d2rivout_in(:,:)=d2rivout_out(:,:)
d2fldout_in(:,:)=d2fldout_out(:,:)
d1pthflw_in(:,:)=d1pthflw_out(:,:)
d2runoff(:,:)=0._jprb

contains




! ==================================================
subroutine storage_sea_surface
use mod_map    ,only: i1next, nseqall
use mod_map    ,only: d2rivlen, d2rivwth, d2elevtn, d2rivelv, d2rivhgt
implicit none

integer(kind=jpim) ::  iseq, jseq, kseq
real(kind=jprb)    ::  dseaelv, dseadph

!!==============================
do iseq=1, nseqall
  jseq=iseq
  do while( i1next(jseq)>0 )
    kseq=jseq
    jseq=i1next(kseq)
  end do
  dseaelv=d2elevtn(jseq,1)

  dseadph=max(dseaelv-d2rivelv(iseq,1),0._jprb)
  dseadph=min(dseadph,d2rivhgt(iseq,1))
  d2rivsto_out(iseq,1)=dseadph*d2rivlen(iseq,1)*d2rivwth(iseq,1)
end do
    
end subroutine storage_sea_surface






! ==================================================
subroutine read_rest_bin
use mod_input ,only: tmpnam
implicit none

real(kind=jprm)     ::  r2temp(nx,ny)
real(kind=jprm)     :: r1pth(npthout,npthlev)
character(len=256)  :: cfile,cdate

integer(kind=jpim)  ::  ipth, ilev, ix, iy

#ifdef UseMPI
character*2         ::  ctmp
#endif
!!==============================

write(cdate,'(i4.4,2i2.2)') iyyyy, imm, idd
cfile=trim(crestdir)//trim(cvnrestart)//trim(cdate)//trim(csufbin)

#ifdef UseMPI
  write(ctmp,'(i2.2)') regionthis
  cfile=cfile//'-'//trim(ctmp)   !
#endif

write(lognam,*)'read_rest: read restart binary: ', trim(cfile)
open(tmpnam,file=cfile,form='unformatted',access='direct',recl=4*nx*ny)
read(tmpnam,rec=1) r2temp
 call map2vec(r2temp,d2rivsto_out)
read(tmpnam,rec=2) r2temp
 call map2vec(r2temp,d2fldsto_out)

if ( .not. lstoonly )then
  read(tmpnam,rec=3) r2temp
   call map2vec(r2temp,d2rivout_pre)
   d2rivout_out=d2rivout_pre
  read(tmpnam,rec=4) r2temp
   call map2vec(r2temp,d2fldout_pre)
   d2fldout_out=d2fldout_pre
  read(tmpnam,rec=5) r2temp
   call map2vec(r2temp,d2rivdph_pre)
  read(tmpnam,rec=6) r2temp
   call map2vec(r2temp,d2fldsto_pre)
endif
close(tmpnam)


if( lpthout )then
  if( .not. lstoonly )then
    cfile=trim(creststo)//'.pth'
    write(lognam,*)'read_rest: read restart binary: ', trim(cfile)

    open(tmpnam,file=cfile,form='unformatted',access='direct',recl=4*npthout*npthlev)
    read(tmpnam,rec=1) r1pth
    d1pthflw_pre(:,:)=r1pth(:,:)
    close(tmpnam)
  endif
endif

end subroutine read_rest_bin





! ==================================================
subroutine read_rest_cdf     !! bugfix v362
#ifdef UseCDF
use netcdf
use lib_netcdf_util
use mod_map ,only: nseqmax
implicit none

integer(kind=jpim)    :: ncid,varid
character(len=256)  :: cfile,cdate
real(kind=jprb)     ::  r2temp(nx,ny)

! yongjiu dai 08/23/2014
!cfile=trim(creststo)
write(cdate,'(i4.4,2i2.2)') iyyyy, imm, idd
cfile=trim(crestdir)//trim(cvnrestart)//trim(cdate)//trim(csufcdf)

write(lognam,*)'read_rest: read restart netcdf: ', trim(cfile)

call ncerror( nf90_open(cfile,nf90_nowrite,ncid), 'opening '//cfile)

call ncerror( nf90_inq_varid(ncid,'rivsto',varid))
call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
call map2vecd(r2temp,d2rivsto_out)

call ncerror( nf90_inq_varid(ncid,'fldsto',varid))
call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
call map2vecd(r2temp,d2fldsto_out)


if( .not. lstoonly )then
  call ncerror( nf90_inq_varid(ncid,'rivout',varid))
  call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
  call map2vecd(r2temp,d2rivout_pre)
  d2rivout_out=d2rivout_pre

  call ncerror( nf90_inq_varid(ncid,'fldout',varid))
  call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
  call map2vecd(r2temp,d2fldout_pre)
  d2fldout_out=d2fldout_pre

  d1pthflw_pre(:,:) = 0._jprb
  d1pthflw_out(:,:)=d1pthflw_pre(:,:)

  call ncerror( nf90_inq_varid(ncid,'rivdph',varid))
  call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
  call map2vecd(r2temp,d2rivdph_pre)

  call ncerror( nf90_inq_varid(ncid,'fldsto_pre',varid))
  call ncerror( nf90_get_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/) ) )
  call map2vecd(r2temp,d2fldsto_pre)

endif

call ncerror( nf90_close(ncid) )

#endif
end subroutine read_rest_cdf








! =================================================
subroutine map2vec(r2temp,d2var)
implicit none

real(kind=jprm)    :: r2temp(nx,ny)
real(kind=jprb)    :: d2var(nseqmax,1)
integer(kind=jpim) :: ix, iy, iseq

do iy=1, ny
  do ix=1, nx
    if( i2vector(ix,iy)>0 )then
      iseq=i2vector(ix,iy)
      d2var(iseq,1)=dble(r2temp(ix,iy))
    endif
  end do
end do

end subroutine map2vec


! =================================================
subroutine map2vecd(r2temp,d2var)  !! bugfix v362
implicit none

real(kind=jprb)    :: r2temp(nx,ny)
real(kind=jprb)    :: d2var(nseqmax,1)
integer(kind=jpim) :: ix, iy, iseq

do iy=1, ny
  do ix=1, nx
    if( i2vector(ix,iy)>0 )then
      iseq=i2vector(ix,iy)
      d2var(iseq,1)=(r2temp(ix,iy))
    endif
  end do
end do

end subroutine map2vecd


end subroutine init_cond
end module init_cond_mod
