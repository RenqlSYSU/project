module control_phy_mod
contains 
subroutine control_phy
! ================================================
! control physics
! ================================================

use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam,dt,dt_def,lkine,ladpstp,lfldout,lpthout,ldam !! **dam
use mod_input  ,only: pgrv, pdstmth, pcadp, irestart,lstoonly
use mod_map    ,only: nseqall,nseqriv,npthout
use mod_prog   ,only: d2rivout_out, d2fldout_out, d1pthflw_out
use mod_diag   ,only: d2rivvel, d2pthout, d2pthinf, d2outflw_avg,d2outflw
use mod_diag   ,only: d2rivout_avg, d2fldout_avg, d2rivvel_avg, d2pthout_avg, d1pthflw_avg, nadd


use calc_fldstg_mod
use calc_rivout_mod
use calc_fldout_mod
use calc_pthout_mod
use calc_damout_mod   !! **dam
use calc_outpre_mod
use calc_stonxt_mod
use calc_watbal_mod
use calc_rivout_kine_mod

implicit none

!! local
integer(kind=jpim)  ::  it, nt
integer(kind=jpim)  ::  ifirst
data                    ifirst / 0 /
save                    ifirst

integer(kind=jpim)            ::  iseq, ipth

!$ save                           iseq, ipth
!!========================================
! write(lognam,*) " "
! write(lognam,*) "****controlphy: ****"

!!=======================================

if( ifirst==0 )then
  write(lognam,*) 'adpstp: daily diagnose'
  call calc_fldstg
  if( lstoonly )then
    call calc_outpre
  endif
  ifirst=1
endif

nt=1
if( ladpstp )then                                                     ! adoptive time step
  call calc_adpstp
endif

!! ==========
do it=1, nt

!!=0 . init

! ====================================
! 1. calculate discharge, diffusive wave
  if ( .not. lkine ) then 
    call calc_rivout

    if( lfldout )then
      call calc_fldout
    endif

    if( lpthout )then
      call calc_pthout
    endif

  else
    call calc_rivout_kine
  endif

  if ( ldam ) then     !! **dam
    call calc_damout
  endif

! ====================================
! 2.  calculate the storage in the next time step in ftcs diff. eq.
  call calc_stonxt

! ====================================
! 3. calculate river and floodplain staging
  call calc_fldstg

! ====================================
! 4.  write water balance monitoring to iofile
  call calc_watbal(it)

! ====================================
! 5. calculate daily average
  nadd=nadd+1
  !$omp parallel do
  do iseq=1, nseqall
    d2rivout_avg(iseq,1)=d2rivout_avg(iseq,1)+d2rivout_out(iseq,1)
    d2fldout_avg(iseq,1)=d2fldout_avg(iseq,1)+d2fldout_out(iseq,1)
    d2rivvel_avg(iseq,1)=d2rivvel_avg(iseq,1)+d2rivvel(iseq,1)
    d2outflw_avg(iseq,1)=d2outflw_avg(iseq,1)+d2outflw(iseq,1)    !! bugfix v362
  end do
  !$omp end parallel do

  if( lpthout )then
    !$omp parallel do
    do iseq=1, nseqall
      d2pthout_avg(iseq,1)=d2pthout_avg(iseq,1)+d2pthout(iseq,1)-d2pthinf(iseq,1)
    end do
    !$omp end parallel do

    !$omp parallel do
    do ipth=1, npthout
      d1pthflw_avg(ipth,:)=d1pthflw_avg(ipth,:)+d1pthflw_out(ipth,:)
    end do
    !$omp end parallel do
  endif

end do

dt=dt_def   !! reset dt



contains

! ==================================================
subroutine calc_adpstp

use mod_input  ,only: lstoonly
use mod_map    ,only: d2nxtdst
use mod_diag   ,only: d2rivdph

implicit none

!$ save
integer(kind=jpim) ::  iseq
real(kind=jprb)    ::  dt_min, ddph, ddst
!$omp threadprivate(ddph,ddst)
! ==========

dt_min=dt_def
!$omp parallel do reduction(min:dt_min)
do iseq=1, nseqriv
  ddph=max(d2rivdph(iseq,1),0.01d0 )
  ddst=d2nxtdst(iseq,1)
  dt_min=min( dt_min, pcadp*ddst * (pgrv*ddph)**(-0.5) )
end do
!$omp end parallel do

!$omp parallel do reduction(min:dt_min)
do iseq=nseqriv+1, nseqall
  ddph=max(d2rivdph(iseq,1),0.01d0 )
  ddst=pdstmth
  dt_min=min( dt_min, pcadp*ddst * (pgrv*ddph)**(-0.5) )
end do
!$omp end parallel do

nt=int( dt_def * dt_min**(-1.) -0.01 )+1
dt=dt_def * real(nt)**(-1.)

if( nt>=2 ) write(lognam,'(a15,i4,3f10.2)') "adpstp: nt=",nt, dt_def, dt_min, dt
! ==================================================
end subroutine calc_adpstp




end subroutine control_phy
end module control_phy_mod
