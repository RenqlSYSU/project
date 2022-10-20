module calc_rivout_mod
contains
subroutine calc_rivout
! ================================================
! calculate discharge, inertial equation
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: dt
use mod_input  ,only: pmanriv,  pdstmth,  pgrv
use mod_map    ,only: i1next,   nseqall,  nseqriv,  nseqmax
use mod_map    ,only: d2rivelv, d2elevtn, d2nxtdst, d2rivwth
use mod_prog   ,only: d2rivsto_out, d2rivout_out
use mod_prog   ,only: d2rivout_pre, d2rivdph_pre
use mod_diag   ,only: d2rivdph, d2rivvel, d2rivinf, d2sfcelv

implicit none

!*** local
      real(kind=jprb)    ::  d2sfcelv_pre(nseqmax,1)                  !! water surface elevation (t-1) [m]
      real(kind=jprb)    ::  d2stoout(nseqmax,1)                      !! total outflow from a grid     [m3]
      real(kind=jprb)    ::  d2rate(nseqmax,1)                        !! outflow correction
!$ save
      integer(kind=jpim) ::  iseq, jseq
      real(kind=jprb)    ::  dslope, darea, dflw, dout_pre, dflw_pre, dflw_imp
!$omp threadprivate   (jseq, dslope, darea, dflw, dout_pre, dflw_pre, dflw_imp)
! ================================================
!$omp parallel do
do iseq=1, nseqall
  d2sfcelv(iseq,1)     = d2rivelv(iseq,1) + d2rivdph(iseq,1)
  d2sfcelv_pre(iseq,1) = d2rivelv(iseq,1) + d2rivdph_pre(iseq,1)
  d2rivinf(iseq,1) = 0.d0
  d2stoout(iseq,1) = 0.d0
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv                                                    !! for normal cells
  jseq=i1next(iseq)

  dslope =  (d2sfcelv(iseq,1)-d2sfcelv(jseq,1)) * d2nxtdst(iseq,1)**(-1.d0)
  dflw = max(d2sfcelv(iseq,1),d2sfcelv(jseq,1)) - d2rivelv(iseq,1)
  darea = d2rivwth(iseq,1) * dflw

  dflw_pre=max(d2sfcelv_pre(iseq,1),d2sfcelv_pre(jseq,1)) - d2rivelv(iseq,1)
  dflw_imp=(dflw*dflw_pre)**0.5d0                                       !! semi implicit flow depth
  if( dflw_imp<=0.d0 ) dflw_imp=dflw

  if( dflw_imp>1.d-5 .and. darea>1.d-5 )then                                                  !! local inertial equation [bates et al., 2010, j.hydro.]
    dout_pre = d2rivout_pre(iseq,1) * d2rivwth(iseq,1)**(-1.d0)                         !! outflow (t-1) [m2/s] (unit width)
    d2rivout_out(iseq,1) = d2rivwth(iseq,1) * ( dout_pre + pgrv*dt*dflw_imp*dslope ) &
                           * ( 1.d0 + pgrv*dt*pmanriv**2.d0*abs(dout_pre)*dflw_imp**(-7.d0/3.d0) )**(-1.d0)
    d2rivvel(iseq,1) = d2rivout_out(iseq,1) * darea**(-1.d0)
  else
    d2rivout_out(iseq,1) = 0.d0
    d2rivvel(iseq,1) = 0.d0
  endif

  if( d2rivout_out(iseq,1) >= 0.d0 )then                                  !! total outflow from a grid
    d2stoout(iseq,1) = d2stoout(iseq,1) + d2rivout_out(iseq,1)*dt
  else
!$omp atomic
    d2stoout(jseq,1) = d2stoout(jseq,1) - d2rivout_out(iseq,1)*dt
  endif
end do
!$omp end parallel do

!$omp parallel do                                                     !! for river mouth grids
do iseq=nseqriv+1, nseqall
  dslope = ( d2sfcelv(iseq,1)-d2elevtn(iseq,1) ) * pdstmth**(-1.d0)
!!  dslope = max( dslope, 0.d0 )                                        !! no backfrom from ocean
  dflw   = d2rivdph(iseq,1)
  darea  = d2rivwth(iseq,1) * dflw

  dflw_pre=d2rivdph_pre(iseq,1)
  dflw_imp=(dflw*dflw_pre)**0.5d0                                     !! semi implicit flow depth
  if( dflw_imp<=0.d0 ) dflw_imp=dflw

  if( dflw_imp>1.d-5 .and. darea>1.d-5 )then 
    dout_pre = d2rivout_pre(iseq,1) * d2rivwth(iseq,1)**(-1.d0)
    d2rivout_out(iseq,1) = d2rivwth(iseq,1) * ( dout_pre + pgrv*dt*dflw_imp*dslope ) &
                           * ( 1.d0 + pgrv*dt*pmanriv**2.d0*abs(dout_pre)*dflw_imp**(-7.d0/3.d0) )**(-1.d0)
    d2rivvel(iseq,1) = d2rivout_out(iseq,1) * darea**(-1.d0)
  else
    d2rivout_out(iseq,1) = 0.d0
    d2rivvel(iseq,1) = 0.d0
  endif
  if( d2rivout_out(iseq,1)>=0.d0 ) d2stoout(iseq,1) = d2stoout(iseq,1) + d2rivout_out(iseq,1)*dt
end do
!$omp end parallel do


!$omp parallel do                                                     !! outflow correcttion if total outflow > storage
do iseq=1, nseqall
  if( d2stoout(iseq,1) > 0.d0)then
    d2rate(iseq,1) = min( d2rivsto_out(iseq,1) * d2stoout(iseq,1)**(-1.d0), 1.d0 )
  else
    d2rate(iseq,1) = 1.d0
  endif
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv
  jseq=i1next(iseq)
  if( d2rivout_out(iseq,1) >= 0.d0 )then
    d2rivout_out(iseq,1) = d2rivout_out(iseq,1)*d2rate(iseq,1)
  else
    d2rivout_out(iseq,1) = d2rivout_out(iseq,1)*d2rate(jseq,1)
  endif
!$omp atomic
  d2rivinf(jseq,1) = d2rivinf(jseq,1) + d2rivout_out(iseq,1)             !! total inflow to a grid (from upstream)

  d2rivout_pre(iseq,1)=d2rivout_out(iseq,1)                              !! save outflow (t)
  d2rivdph_pre(iseq,1)=d2rivdph(iseq,1)                              !! save depth   (t)
end do
!$omp end parallel do

!$omp parallel do
do iseq=nseqriv+1, nseqall
  d2rivout_out(iseq,1) = d2rivout_out(iseq,1)*d2rate(iseq,1)

  d2rivout_pre(iseq,1)=d2rivout_out(iseq,1)                              !! save outflow (t)
  d2rivdph_pre(iseq,1)=d2rivdph(iseq,1)                              !! save depth   (t)
end do
!$omp end parallel do

end subroutine calc_rivout
end module calc_rivout_mod
