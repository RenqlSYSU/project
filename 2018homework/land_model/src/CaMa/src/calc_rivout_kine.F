module calc_rivout_kine_mod
contains
subroutine calc_rivout_kine
! ================================================
! calculate discharge, kinematic wave
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: dt
use mod_input  ,only: pmanriv, pdstmth, pminslp
use mod_map    ,only: i1next,   nseqall,  nseqriv
use mod_map    ,only: d2rivelv, d2elevtn, d2nxtdst, d2rivwth
use mod_prog   ,only: d2rivsto_out, d2rivout_out
use mod_diag   ,only: d2rivdph, d2rivvel, d2rivinf, d2sfcelv

implicit none

!$ save
      integer(kind=jpim) ::  iseq, jseq
      real(kind=jprb)    ::  dslope, darea, dvel
!$omp threadprivate   (jseq, dslope, darea, dvel)
!
! ================================================
!$omp parallel do
do iseq=1, nseqall
  d2sfcelv(iseq,1) = d2rivelv(iseq,1) + d2rivdph(iseq,1)          !! diffusive
  d2rivinf(iseq,1) = 0.d0
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv
  jseq   = i1next(iseq)
  dslope = (d2elevtn(iseq,1)-d2elevtn(jseq,1)) * d2nxtdst(iseq,1)**(-1.)
  dslope = max(dslope,pminslp)
  dvel   = pmanriv**(-1.) * dslope**0.5 * d2rivdph(iseq,1)**(2./3.)
  darea  = d2rivwth(iseq,1) * d2rivdph(iseq,1)

  d2rivvel(iseq,1) = dvel
  d2rivout_out(iseq,1) = darea * dvel
  d2rivout_out(iseq,1) = min(  d2rivout_out(iseq,1), d2rivsto_out(iseq,1)/dt )

  !$omp atomic
  d2rivinf(jseq,1) = d2rivinf(jseq,1) + d2rivout_out(iseq,1)
end do
!$omp end parallel do

!$omp parallel do
do iseq=nseqriv+1, nseqall
  dslope = pminslp
  dvel   = pmanriv**(-1.) * dslope**0.5 * d2rivdph(iseq,1)**(2./3.)
  darea  = d2rivwth(iseq,1) * d2rivdph(iseq,1)

  d2rivvel(iseq,1) = dvel
  d2rivout_out(iseq,1) = darea * dvel
  d2rivout_out(iseq,1) = min(  d2rivout_out(iseq,1), d2rivsto_out(iseq,1)/dt )
end do
!$omp end parallel do


end subroutine calc_rivout_kine
end module calc_rivout_kine_mod
