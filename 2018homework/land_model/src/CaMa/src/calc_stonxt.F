module calc_stonxt_mod
contains
subroutine calc_stonxt
! ================================================
! calculate the storage in the next time step in ftcs diff. eq.
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: dt
use mod_map    ,only: nseqall
use mod_prog   ,only: d2rivout_out,d2fldout_out,d2rivsto_out,d2fldsto_out,d2runoff
use mod_diag   ,only: d2rivinf,    d2fldinf,    d2pthout,    d2pthinf,    d2fldfrc,    d2outflw,    d2storge
use mod_diag   ,only: dglbstopre,  dglbstonxt,  dglbstonew,  dglbrivinf,  dglbrivout
implicit none

!*** local
!$ save
integer(kind=jpim)    ::  iseq
real(kind=jprb)       ::  drivrof, dfldrof
!$omp threadprivate      (drivrof, dfldrof)
!!==============================
dglbstopre=0.d0
dglbstonxt=0.d0
dglbstonew=0.d0
dglbrivinf=0.d0
dglbrivout=0.d0

!$omp parallel do reduction(+:dglbstopre,dglbrivinf,dglbrivout,dglbstonxt,dglbstonew)
do iseq=1, nseqall

  dglbstopre = dglbstopre + d2rivsto_out(iseq,1)    + d2fldsto_out(iseq,1)
  dglbrivinf = dglbrivinf + d2rivinf(iseq,1)*dt     + d2fldinf(iseq,1)*dt     + d2pthinf(iseq,1)*dt
  dglbrivout = dglbrivout + d2rivout_out(iseq,1)*dt + d2fldout_out(iseq,1)*dt + d2pthout(iseq,1)*dt

  d2rivsto_out(iseq,1) = d2rivsto_out(iseq,1) + d2rivinf(iseq,1)*dt - d2rivout_out(iseq,1)*dt
  d2rivsto_out(iseq,1) = max( d2rivsto_out(iseq,1), 0.d0 )

  d2fldsto_out(iseq,1) = d2fldsto_out(iseq,1) + d2fldinf(iseq,1)*dt - d2fldout_out(iseq,1)*dt &
                                              + d2pthinf(iseq,1)*dt - d2pthout(iseq,1)*dt
  if( d2fldsto_out(iseq,1) < 0.d0 )then
    d2rivsto_out(iseq,1)=max( d2rivsto_out(iseq,1)+d2fldsto_out(iseq,1), 0.d0 )
    d2fldsto_out(iseq,1)=0.d0
  endif

  dglbstonxt = dglbstonxt + d2rivsto_out(iseq,1) + d2fldsto_out(iseq,1)
  d2outflw(iseq,1)=d2rivout_out(iseq,1)+d2fldout_out(iseq,1)+d2pthout(iseq,1)

  drivrof = d2runoff(iseq,1) * (1.d0-d2fldfrc(iseq,1)) * dt
  dfldrof = d2runoff(iseq,1) *       d2fldfrc(iseq,1)  * dt
  d2rivsto_out(iseq,1) = d2rivsto_out(iseq,1) + drivrof
  d2fldsto_out(iseq,1) = d2fldsto_out(iseq,1) + dfldrof

  dglbstonew=dglbstonew+d2rivsto_out(iseq,1)+d2fldsto_out(iseq,1)
  d2storge(iseq,1)=d2rivsto_out(iseq,1)+d2fldsto_out(iseq,1)

end do
!$omp end parallel do

end subroutine calc_stonxt
end module calc_stonxt_mod
