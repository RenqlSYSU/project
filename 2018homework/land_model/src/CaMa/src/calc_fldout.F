module calc_fldout_mod
contains
subroutine calc_fldout
! ================================================
! calculation of floodplain flow
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: dt
use mod_input  ,only: pmanfld,  pdstmth,  pgrv
use mod_map    ,only: i1next,   nseqall,  nseqriv,  nseqmax
use mod_map    ,only: d2rivelv, d2elevtn, d2nxtdst, d2rivlen, d2rivhgt, d2rivwth
use mod_prog   ,only: d2fldsto_out, d2fldout_out, d2rivout_out
use mod_prog   ,only: d2fldout_pre, d2rivdph_pre, d2fldsto_pre
use mod_diag   ,only: d2flddph, d2fldinf, d2sfcelv

implicit none

!*** local
      real(kind=jprb)    ::  d2sfcelv_now(nseqmax,1)                  !! water surface elev (t)  [m] (for floodplain)
      real(kind=jprb)    ::  d2sfcelv_pre(nseqmax,1)                  !! water surface elev (t-1)[m] (for floodplain)
      real(kind=jprb)    ::  d2flware_now(nseqmax,1)                  !! flow area (t)           [m2]
      real(kind=jprb)    ::  d2flware_pre(nseqmax,1)                  !! flow area (t-1)         [m2]
      real(kind=jprb)    ::  d2flddph_pre(nseqmax,1)                  !! floodplain water depth  [m]
      real(kind=jprb)    ::  d2stoout(nseqmax,1)                      !! total outflow from grid [m3]
      real(kind=jprb)    ::  d2rate(nseqmax,1)                        !! outflow correction
!$ save
      real(kind=jprb)    ::  dslope, dflw, dout_pre, dflw_pre, dflw_imp, dare_imp
      integer(kind=jpim) ::  iseq, jseq
!$omp threadprivate         (dslope, dflw, dout_pre, dflw_pre, dflw_imp, dare_imp, jseq)
! ================================================
!$omp parallel do
do iseq=1, nseqall
  d2sfcelv_now(iseq,1) = max( d2elevtn(iseq,1), d2sfcelv(iseq,1) )
  d2sfcelv_pre(iseq,1) = max( d2elevtn(iseq,1), d2rivelv(iseq,1)+d2rivdph_pre(iseq,1) )
  d2flddph_pre(iseq,1) = max( d2rivdph_pre(iseq,1)-d2rivhgt(iseq,1), 0.d0 )
  d2flware_now(iseq,1) = d2fldsto_out(iseq,1) * d2rivlen(iseq,1)**(-1.d0)
  d2flware_pre(iseq,1) = d2fldsto_pre(iseq,1) * d2rivlen(iseq,1)**(-1.d0)
  d2flware_now(iseq,1) = max( d2flware_now(iseq,1) - d2flddph(iseq,1)    *d2rivwth(iseq,1), 0.d0 ) !! remove river channel area
  d2flware_pre(iseq,1) = max( d2flware_pre(iseq,1) - d2flddph_pre(iseq,1)*d2rivwth(iseq,1), 0.d0 )
  d2fldinf(iseq,1) = 0.d0
  d2stoout(iseq,1) = 0.d0
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv                                                    !! for normal cell
  jseq=i1next(iseq)

  dslope  = (d2sfcelv_now(iseq,1)-d2sfcelv_now(jseq,1)) * d2nxtdst(iseq,1)**(-1.d0)

  dflw = min(d2sfcelv_now(iseq,1),d2sfcelv_now(jseq,1)) - d2elevtn(iseq,1)        !! should be min[i,j] to avoid instability
  dflw = max(dflw,0.d0)
  dflw_pre = min(d2sfcelv_pre(iseq,1),d2sfcelv_pre(jseq,1)) - d2elevtn(iseq,1)
  dflw_pre = max(dflw_pre,0.d0)

  dflw_imp = (dflw*dflw_pre)**0.5d0                                     !! semi implicit flow depth

  dare_imp =      (d2flware_now(iseq,1)*d2flware_pre(iseq,1))**0.5d0    !! semi implicit flow area
  dare_imp = min( (d2flware_now(jseq,1)*d2flware_pre(jseq,1))**0.5d0, dare_imp )

  if( dflw_imp<=0.d0 ) dflw_imp=dflw
  if( dare_imp<=0.d0 ) dare_imp=d2flware_now(iseq,1)

  if( dflw_imp>1.d-5 .and. dare_imp>1.d-5 )then                         !! local inertial equation, see [bates et al., 2010, j.hydrol.]
    dout_pre = d2fldout_pre(iseq,1)
    d2fldout_out(iseq,1) = ( dout_pre + pgrv*dt*dare_imp*dslope ) &
                           * (1.d0 + pgrv*dt*pmanfld**2.d0*abs(dout_pre)*dflw_imp**(-4.d0/3.d0)*dare_imp**(-1.d0) )**(-1.d0)
  else
    d2fldout_out(iseq,1) = 0.d0
  endif
  if( d2fldout_out(iseq,1)*d2rivout_out(iseq,1)<0.d0 ) d2fldout_out(iseq,1)=0.d0

  if( d2fldout_out(iseq,1) >= 0.d0 )then                                  !! total outflow from each grid
    d2stoout(iseq,1) = d2stoout(iseq,1) + d2fldout_out(iseq,1)*dt
  else
!$omp atomic
    d2stoout(jseq,1) = d2stoout(jseq,1) - d2fldout_out(iseq,1)*dt
  endif
end do
!$omp end parallel do

!$omp parallel do
do iseq=nseqriv+1, nseqall                                            !! for river mouth cell
  dslope = ( d2sfcelv_now(iseq,1)-d2elevtn(iseq,1) ) * pdstmth**(-1.d0)
!  dslope = max( dslope, 0.d0 )                                        !! no backward flow from ocean
  dflw   = d2flddph(iseq,1)

  dflw_pre=max(d2sfcelv_pre(iseq,1)-d2elevtn(iseq,1),0.d0)
  dflw_imp=(dflw*dflw_pre)**0.5d0
  dare_imp = (d2flware_now(iseq,1)*d2flware_pre(iseq,1))**0.5d0
  if( dflw_imp<=0.d0 ) dflw_imp=dflw
  if( dare_imp<=0.d0 ) dare_imp=d2flware_now(iseq,1)

  if( dflw_imp>1.d-5 .and. dare_imp>1.d-5 )then
    dout_pre = d2fldout_pre(iseq,1)
    d2fldout_out(iseq,1) = ( dout_pre + pgrv*dt*dare_imp*dslope ) &
                           * (1.d0 + pgrv*dt*pmanfld**2.d0*abs(dout_pre)*dflw_imp**(-4.d0/3.d0)*dare_imp**(-1.d0) )**(-1.d0)
  else
    d2fldout_out(iseq,1) = 0.d0
  endif
  if( d2fldout_out(iseq,1) >= 0.d0 ) d2stoout(iseq,1) = d2stoout(iseq,1) + d2fldout_out(iseq,1)*dt
end do
!$omp end parallel do


!$omp parallel do                                                     !! calculate total outflow from a grid
do iseq=1, nseqall                                                    !! outflow correction applied when the total outflow > storage
  if( d2stoout(iseq,1) > 0.d0 )then
    d2rate(iseq,1) = min( d2fldsto_out(iseq,1) * d2stoout(iseq,1)**(-1.), 1.d0 )
  else
    d2rate(iseq,1) = 1.d0
  endif
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv
  jseq=i1next(iseq)
  if( d2fldout_out(iseq,1) >= 0.d0 )then
    d2fldout_out(iseq,1) = d2fldout_out(iseq,1)*d2rate(iseq,1)
  else
    d2fldout_out(iseq,1) = d2fldout_out(iseq,1)*d2rate(jseq,1)
  endif
!$omp atomic
  d2fldinf(jseq,1) = d2fldinf(jseq,1) + d2fldout_out(iseq,1)             !! total inflow [m3/s] (from upstream)

  d2fldout_pre(iseq,1)=d2fldout_out(iseq,1)
  d2fldsto_pre(iseq,1)=d2fldsto_out(iseq,1)
end do
!$omp end parallel do

!$omp parallel do
do iseq=nseqriv+1, nseqall
  d2fldout_out(iseq,1) = d2fldout_out(iseq,1)*d2rate(iseq,1)

  d2fldout_pre(iseq,1)=d2fldout_out(iseq,1)
  d2fldsto_pre(iseq,1)=d2fldsto_out(iseq,1)
end do
!$omp end parallel do

end subroutine calc_fldout
end module calc_fldout_mod
