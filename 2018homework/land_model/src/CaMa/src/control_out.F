module control_out_mod
contains 
subroutine control_out(iyyyyp,immp,iddp)

! ================================================
! control output files
! ================================================
use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam, nx,ny, imis,rmis, &
                      loutcdf, loutvec, lpthout, coutinsdir
use mod_map    ,only: nseqall, nseqmax, npthout, npthlev
use mod_output ,only: funits,  cdfunits, coutvardir, coutvn, nvars
use mod_time   ,only: kstep, kyyyymmdd, khhmm, khhmmp
use mod_prog   ,only: d2rivsto_out, d2fldsto_out
use mod_diag   ,only: d2rivdph, d2flddph, d2fldfrc, d2fldare,&
                      d2sfcelv, d2outflw, d2storge, d2outflw_avg
use mod_diag   ,only: d2rivout_avg, d2fldout_avg, d2rivvel_avg, d2pthout_avg, d1pthflw_avg, nadd, d2outins

use create_outbin_mod
use create_outcdf_mod

implicit none
! ================================================
!yongjiu dai 08/23/2014
!integer(kind=jpim)    ::    iyyyyp,immp,iddp     !! date (start of time step)
integer(kind=jpim),intent(in)  :: iyyyyp          !! date
integer(kind=jpim),intent(in)  :: immp,iddp       !! date (previous)

save
integer(kind=jpim)             :: irec
integer(kind=jpim)             :: iseq, ipth
! ================================================


if ( kstep == 1 .or. (immp==1 .and. iddp==1 .and. khhmmp==0) ) then
  if (.not.loutcdf)then
    call create_outbin(iyyyyp)
  else
    call create_outcdf(iyyyyp)
  endif
  irec=1
endif

if (khhmm==0 ) then
  write(lognam,*),'control out write',nadd,irec
  !$omp parallel do
  do iseq=1, nseqall
    d2rivout_avg(iseq,1)=d2rivout_avg(iseq,1)/dble(nadd)
    d2fldout_avg(iseq,1)=d2fldout_avg(iseq,1)/dble(nadd)
    d2rivvel_avg(iseq,1)=d2rivvel_avg(iseq,1)/dble(nadd)
    d2outflw_avg(iseq,1)=d2outflw_avg(iseq,1)/dble(nadd)   !! bugfix v362
  end do
  !$omp end parallel do

  if( lpthout )then
    !$omp parallel do
    do iseq=1, nseqall
      d2pthout_avg(iseq,1)=d2pthout_avg(iseq,1)/dble(nadd)
    end do
    !$omp end parallel do

    !$omp parallel do
    do ipth=1, npthout
      d1pthflw_avg(ipth,:)=d1pthflw_avg(ipth,:)/dble(nadd)
    end do
    !$omp end parallel do
  endif
  nadd=0

  if( coutinsdir .ne. 'NONE' ) call calc_outins
  if ( .not.loutcdf ) then
    call wrte_bin
  else
    call wrte_cdf
  endif 
  irec=irec+1
  if( mod(kyyyymmdd,10000)==101 .and. khhmm==0 ) irec=1

  !$omp parallel do
  do iseq=1, nseqall
    d2rivout_avg(iseq,1)=0.d0
    d2fldout_avg(iseq,1)=0.d0
    d2rivvel_avg(iseq,1)=0.d0
    d2outflw_avg(iseq,1)=0.d0    !! bugfix v362
    nadd=0
  end do
  !$omp end parallel do

  if( lpthout )then
    !$omp parallel do
    do iseq=1, nseqall
      d2pthout_avg(iseq,1)=0.d0
    end do
    !$omp end parallel do

    !$omp parallel do
    do ipth=1, npthout
      d1pthflw_avg(ipth,:)=0.d0
    end do
    !$omp end parallel do
  endif

endif 

contains 





subroutine wrte_bin
! =================================================
implicit none

if (.not.loutvec) then ! vectorial output
  call wrte_outbin(funits(1), irec,d2rivout_avg,coutvardir(1))
  call wrte_outbin(funits(2), irec,d2rivsto_out,coutvardir(2))
  call wrte_outbin(funits(3), irec,d2rivdph,    coutvardir(3))
  call wrte_outbin(funits(4), irec,d2rivvel_avg,coutvardir(4))

  call wrte_outbin(funits(5), irec,d2fldout_avg,coutvardir(5))
  call wrte_outbin(funits(6), irec,d2fldsto_out,coutvardir(6))
  call wrte_outbin(funits(7), irec,d2flddph,    coutvardir(7))
  call wrte_outbin(funits(8), irec,d2fldfrc,    coutvardir(8))
  call wrte_outbin(funits(9), irec,d2fldare,    coutvardir(9))

  call wrte_outbin(funits(10),irec,d2sfcelv,    coutvardir(10))
  call wrte_outbin(funits(11),irec,d2outflw_avg,coutvardir(11))    !! bugfix v362
  call wrte_outbin(funits(12),irec,d2storge,    coutvardir(12))
  call wrte_outbin(funits(13),irec,d2storge,    coutvardir(13))

  call wrte_outbin(funits(14),irec,d2pthout_avg,coutvardir(14))
  call wrte_outpth(funits(15),irec,d1pthflw_avg,coutvardir(15))

else
  call wrte_outvec(funits(1), irec,d2rivout_avg,coutvardir(1))
  call wrte_outvec(funits(2), irec,d2rivsto_out,coutvardir(2))
  call wrte_outvec(funits(3), irec,d2rivdph,    coutvardir(3))
  call wrte_outvec(funits(4), irec,d2rivvel_avg,coutvardir(4))

  call wrte_outvec(funits(5), irec,d2fldout_avg,coutvardir(5))
  call wrte_outvec(funits(6), irec,d2fldsto_out,coutvardir(6))
  call wrte_outvec(funits(7), irec,d2flddph,    coutvardir(7))
  call wrte_outvec(funits(8), irec,d2fldfrc,    coutvardir(8))
  call wrte_outvec(funits(9), irec,d2fldare,    coutvardir(9))

  call wrte_outvec(funits(10),irec,d2sfcelv,    coutvardir(10))
  call wrte_outvec(funits(11),irec,d2outflw_avg,coutvardir(11))    !! bugfix v362
  call wrte_outvec(funits(12),irec,d2storge,    coutvardir(12))
  call wrte_outvec(funits(13),irec,d2outins,    coutvardir(13))

  call wrte_outvec(funits(14),irec,d2pthout_avg,coutvardir(14))
  call wrte_outpth(funits(15),irec,d1pthflw_avg,coutvardir(15))
endif

end subroutine wrte_bin





! ================================================
subroutine wrte_outbin(ifn,irec,d2vec,cvaroutdir)
implicit none

!*** input
integer(kind=jpim),intent(in)  ::  ifn                 !! file number
integer(kind=jpim),intent(in)  ::  irec                !! record
real(kind=jprb),intent(in)     ::  d2vec(nseqmax,1)    !! output data
character*128                  ::  cvaroutdir          !! variable output dir

!*** local
real(kind=jprm)                ::  r2map(nx,ny)
!!==============================
if (cvaroutdir.ne.'NONE') then
  call vec2map(d2vec,r2map)
  write(ifn,rec=irec) r2map
endif
!
end subroutine wrte_outbin




! ================================================
subroutine wrte_outvec(ifn,irec,d2outdat,cvaroutdir)
implicit none

!*** input
integer(kind=jpim),intent(in)             ::  ifn                 !! file number
integer(kind=jpim),intent(in)             ::  irec                !! record
real(kind=jprb),intent(in)                ::  d2outdat(nseqmax,1) !! output data
character*128       ::  cvaroutdir                                !! variable output dir

!*** local
real(kind=jprm)     ::  r2outdat(nseqmax,1)
! ================================================
if (cvaroutdir.ne.'NONE') then
  r2outdat(:,:)=real(d2outdat(:,:))
  write(ifn,rec=irec) r2outdat
endif

end subroutine wrte_outvec



! ================================================
subroutine wrte_outpth(ifn,irec,d2outdat,cvaroutdir)
implicit none

!*** input
integer(kind=jpim),intent(in)             ::  ifn                 !! file number
integer(kind=jpim),intent(in)             ::  irec                !! record
real(kind=jprb),intent(in)                ::  d2outdat(npthout,npthlev) !! output data
character*128       ::  cvaroutdir                                !! variable output dir

!*** local
real(kind=jprm)     ::  r2outdat(npthout,npthlev)
! ================================================
if (cvaroutdir.ne.'NONE') then
  r2outdat(:,:)=real(d2outdat(:,:))
  write(ifn,rec=irec) r2outdat
endif

end subroutine wrte_outpth




subroutine wrte_cdf

#ifdef UseCDF
! ==================================================
use netcdf 
use lib_netcdf_util
use mod_time, only:kminstart,kmin,nsteps
implicit none

integer(kind=jpim)   :: jf,varid
real(kind=jprb)      :: xtime ! seconds since start of the run !
real(kind=jprb)      :: d2vec(nseqmax,1)
real(kind=jprm)      :: z2map(nx,ny)

xtime=real( (kmin-kminstart),jprb) *60._jprb

do jf=1,nvars-1
  if ( cdfunits(jf) .eq. 0 ) cycle 

  select case (jf)
    case (1);  d2vec(:,:)=d2rivout_avg(:,:)
    case (2);  d2vec(:,:)=d2rivsto_out(:,:)
    case (3);  d2vec(:,:)=d2rivdph(:,:)
    case (4);  d2vec(:,:)=d2rivvel_avg(:,:)

    case (5);  d2vec(:,:)=d2fldout_avg(:,:)
    case (6);  d2vec(:,:)=d2fldsto_out(:,:)
    case (7);  d2vec(:,:)=d2flddph(:,:)
    case (8);  d2vec(:,:)=d2fldfrc(:,:)
    case (9);  d2vec(:,:)=d2fldare(:,:)

    case (10); d2vec(:,:)=d2sfcelv(:,:)
    case (11); d2vec(:,:)=d2outflw_avg(:,:)  !! bugfix v362
    case (12); d2vec(:,:)=d2storge(:,:)
    case (13); d2vec(:,:)=d2outins(:,:)
    
    case (14); d2vec(:,:)=d2pthout_avg(:,:)
  end select

  call ncerror( nf90_inq_varid(cdfunits(jf),'time',varid))
  call ncerror( nf90_put_var(cdfunits(jf),varid,xtime,(/irec/)) )

  call ncerror( nf90_inq_varid(cdfunits(jf),trim(coutvn(jf)),varid))

  if (.not. loutvec) then
    call vec2map(d2vec,z2map)
    print *,"wrte_cdf",trim(coutvn(jf)),(/1,1,irec/),(/nx,ny,1/)
    call ncerror( nf90_put_var(cdfunits(jf),varid,z2map(1:nx,1:ny),(/1,1,irec/),(/nx,ny,1/)) )
  else
    print *,"wrte_cdf",trim(coutvn(jf)),(/1,irec/),(/nseqmax,1/)
    call ncerror( nf90_put_var(cdfunits(jf),varid,d2vec(1:nseqmax,1),(/1,irec/),(/nseqmax,1/)) )
  endif
    
  call ncerror( nf90_sync(cdfunits(jf)) )
  
  if ( kstep == nsteps) call ncerror( nf90_close(cdfunits(jf)) )
end do
#endif
end subroutine wrte_cdf 





! ================================================
subroutine vec2map(d2vec,r2map)
use mod_map    ,only: i2vector
implicit none

real(kind=jprb)              :: d2vec(nseqmax,1)
real(kind=jprm)              :: r2map(nx,ny)

!$ save
integer(kind=jpim)           ::  ix,iy,iseq
!$omp threadprivate(ix,iseq)
! ======
r2map(:,:) = rmis
!$omp parallel do
do iy=1, ny
  do ix=1, nx
    if (i2vector(ix,iy)>0) then
      iseq=i2vector(ix,iy)
      r2map(ix,iy) = real(d2vec(iseq,1))
    end if
  end do
end do
!$omp end parallel do

end subroutine vec2map





! ================================================
subroutine calc_outins
! to   calculate discharge, inst. no routing delay

use mod_map    ,only: i2rivseq,rivseqmax,i2vector,i2nextx,i2nexty
use mod_prog   ,only: d2runoff
implicit none

!!** local
integer(kind=jpim)          :: irivseq,iseq,jseq,ix,iy,jx,jy

!!==============================
d2outins(:,:)=d2runoff(:,:)

do irivseq=1, rivseqmax
  do iy=1, ny
    do ix=1, nx
      if( i2rivseq(ix,iy)==irivseq .and. i2nextx(ix,iy)>0 )then
        jx=i2nextx(ix,iy)
        jy=i2nexty(ix,iy)
        iseq=i2vector(ix,iy)
        jseq=i2vector(jx,jy)
        d2outins(jseq,1)=d2outins(jseq,1)+d2outins(iseq,1)
      endif
    end do
  end do
end do

end subroutine calc_outins



end subroutine control_out
end module control_out_mod
