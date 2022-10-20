module control_inp_mod

! modifed by yongjiu dai for coupling with CoLM

contains 

subroutine control_inp(r2roffin)

! =======================================================
! control input runoff
! =======================================================

use parkind1  ,only: jpim, jprb, jprm
use mod_input ,only: nx, ny, nxin, nyin, inpn, linterp, rmis
use mod_map   ,only: i1seqx, i1seqy, nseqall, d2grarea
use mod_prog  ,only: d2runoff

IMPLICIT NONE

real(kind=jprm), intent(in) :: r2roffin(nxin,nyin)  ! for input runoff (mm/day)

!$ save
integer(kind=jpim)            ::  iseq

! ===============================================
!$omp parallel do
do iseq=1, nseqall
  d2runoff(iseq,1)=0.d0
enddo
!$omp end parallel do

! interpolation & unit conversion (mm/day --> m3/sec)
if( linterp )then
  call roff_interp(r2roffin,d2runoff)
else
  call conv_resol(r2roffin,d2runoff)
endif



contains



! ================================================
subroutine roff_interp(r2roffin,d2runoff)
use mod_map    ,only: i1seqx, i1seqy, nseqall
use mod_map    ,only: inpx, inpy, inpa

IMPLICIT NONE

real(kind=jprm),intent(in)  ::  r2roffin(:,:)     !! mm/day
real(kind=jprb),intent(out) ::  d2runoff(:,:)     !! m3/s

!$ save
integer(kind=jpim)  ::  iseq
integer(kind=jpim)  ::  ix, iy, ixin, iyin, inpi  !! for output
!$omp threadprivate    (ix, iy, ixin, iyin, inpi)
! ================================================
!$omp parallel do
do iseq=1, nseqall
  ix=i1seqx(iseq)
  iy=i1seqy(iseq)
  d2runoff(iseq,1)=0.d0
  do inpi=1, inpn
    ixin=inpx(ix,iy,inpi)
    iyin=inpy(ix,iy,inpi)
    if( ixin <= nxin .and. iyin <= nyin )then  ! /04/2014/ yongjiu dai
      if( ixin>0 )then
        if( r2roffin(ixin,iyin).ne.rmis )then
          d2runoff(iseq,1) = d2runoff(iseq,1) + r2roffin(ixin,iyin) * inpa(ix,iy,inpi) * 1.d-3 / 86400.
          d2runoff(iseq,1) = max(d2runoff(iseq,1), 0.d0)
        endif
      endif
    endif
  end do
  d2runoff(iseq,1)=max(d2runoff(iseq,1), 0.d0)
end do
!$omp end parallel do

return
end subroutine roff_interp



! ================================================
subroutine conv_resol(r2roffin,d2runoff)
use mod_map    ,only: i1seqx, i1seqy, nseqall
use mod_map    ,only: d2grarea

IMPLICIT NONE

real(kind=jprm),intent(in)  ::  r2roffin(:,:)     !! mm/day
real(kind=jprb),intent(out) ::  d2runoff(:,:)     !! m3/s

integer(kind=jpim)  ::  by

!$ save
integer(kind=jpim)  ::  iseq
integer(kind=jpim)  ::  ixin, iyin
!$omp threadprivate    (ixin, iyin)
! ================================================
by = nx/nxin
!$omp parallel do
do iseq=1, nseqall
  ixin=int((i1seqx(iseq)-1)/by)+1
  iyin=int((i1seqy(iseq)-1)/by)+1
  if( ixin <= nxin .and. iyin <= nyin )then  ! /04/2014/ yongjiu dai
    if( r2roffin(ixin,iyin).ne.rmis )then
      d2runoff(iseq,1) = r2roffin(ixin,iyin) * d2grarea(iseq,1) * 1.d-3 / 86400.
      d2runoff(iseq,1) = max(d2runoff(iseq,1), 0.d0)
    else
      d2runoff(iseq,1)=0.d0
    endif
  endif
end do
!$omp end parallel do

return
end subroutine conv_resol


end subroutine control_inp
end module control_inp_mod
! --------- EOP ----------
