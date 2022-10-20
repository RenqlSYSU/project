module calc_pthout_mod
contains
subroutine calc_pthout
! ================================================
! calculation of floodplain flow
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: dt
use mod_input  ,only: pgrv
use mod_map    ,only: nseqmax, npthout, npthlev, pth_upst, pth_down, pth_dst, pth_elv, pth_wth, pth_man
use mod_map    ,only: d2rivelv
use mod_prog   ,only: d2rivsto_out, d2fldsto_out, d1pthflw_out, d2rivout_out, d2fldout_out
use mod_prog   ,only: d1pthflw_pre, d2rivdph_pre
use mod_diag   ,only: d2pthout, d2pthinf, d2rivinf, d2fldinf, d2sfcelv

implicit none

!*** local
      real(kind=jprb)    ::  d2sfcelv_pre(nseqmax,1)                  !! water surface elev (t-1)[m] (for floodplain)
      real(kind=jprb)    ::  d2rate(nseqmax,1)                        !! outflow correction

!$ save
      integer(kind=jpim) ::  ipth, ilev, iseq, iseqp, jseqp
      real(kind=jprb)    ::  dslope, dflw, dout_pre, dflw_pre, dflw_imp, dsto_tmp
!$omp threadprivate         (dslope, dflw, dout_pre, dflw_pre, dflw_imp, dsto_tmp, ilev, iseqp, jseqp)
! ================================================
!$omp parallel do
do iseq=1, nseqmax
  d2sfcelv_pre(iseq,1) = d2rivelv(iseq,1)+d2rivdph_pre(iseq,1)
  d2pthout(iseq,1) = 0.d0
  d2pthinf(iseq,1) = 0.d0
  d2rate(iseq,1)=-999.
end do
!$omp end parallel do

!$omp parallel do
do ipth=1, npthout  
  iseqp=pth_upst(ipth)
  jseqp=pth_down(ipth)

  dslope  = (d2sfcelv(iseqp,1)-d2sfcelv(jseqp,1)) * pth_dst(ipth)**(-1.d0)

  do ilev=1, npthlev

    dflw = max(d2sfcelv(iseqp,1),d2sfcelv(jseqp,1)) - pth_elv(ipth,ilev) 
    dflw = max(dflw,0.d0)

    dflw_pre = min(d2sfcelv_pre(iseqp,1),d2sfcelv_pre(jseqp,1)) - pth_elv(ipth,ilev)
    dflw_pre = max(dflw_pre,0.d0)
    dflw_imp = (dflw*dflw_pre)**0.5d0                                       !! semi implicit flow depth

    if( dflw_imp<=0.d0 ) dflw_imp=dflw

    if( dflw_imp>1.d-5 )then                         !! local inertial equation, see [bates et al., 2010, j.hydrol.]
      dout_pre = d1pthflw_pre(ipth,ilev) * pth_wth(ipth,ilev)**(-1.d0)                         !! outflow (t-1) [m2/s] (unit width)
      d1pthflw_out(ipth,ilev) = pth_wth(ipth,ilev) * ( dout_pre + pgrv*dt*dflw_imp*dslope ) &
                         * ( 1.d0 + pgrv*dt*pth_man(ilev)**2.d0*abs(dout_pre)*dflw_imp**(-7.d0/3.d0) )**(-1.d0)
    else
      d1pthflw_out(ipth,ilev) = 0.d0
    endif

    if( d1pthflw_out(ipth,ilev) >= 0.d0 )then                                  !! total outflow from each grid
!$omp atomic
      d2pthout(iseqp,1) = d2pthout(iseqp,1) + d1pthflw_out(ipth,ilev)
    else
!$omp atomic
      d2pthout(jseqp,1) = d2pthout(jseqp,1) - d1pthflw_out(ipth,ilev)
    endif
  end do
end do
!$omp end parallel do

!$omp parallel do                                              !! calculate total outflow from a grid
do iseq=1, nseqmax
  if( d2pthout(iseq,1) > 1.d-10 )then
    dsto_tmp = ( d2rivsto_out(iseq,1)+d2fldsto_out(iseq,1) ) &
                  - d2rivout_out(iseq,1)*dt + d2rivinf(iseq,1)*dt - d2fldout_out(iseq,1)*dt + d2fldinf(iseq,1)*dt
    d2rate(iseq,1) = min( dsto_tmp * (d2pthout(iseq,1)*dt)**(-1.), 1.d0 )
  else
    d2rate(iseq,1) = 1.d0
  endif
  d2pthout(iseq,1) = d2pthout(iseq,1) * d2rate(iseq,1)
end do
!$omp end parallel do

!$omp parallel do   
do ipth=1, npthout
  iseqp=pth_upst(ipth)
  jseqp=pth_down(ipth)

  do ilev=1, npthlev
    if( d1pthflw_out(ipth,ilev) >= 0.d0 )then
      d1pthflw_out(ipth,ilev) = d1pthflw_out(ipth,ilev)*d2rate(iseqp,1)
!$omp atomic
      d2pthinf(jseqp,1) = d2pthinf(jseqp,1) + d1pthflw_out(ipth,ilev)             !! total inflow [m3/s] (from upstream)
    else
      d1pthflw_out(ipth,ilev) = d1pthflw_out(ipth,ilev)*d2rate(jseqp,1)
!$omp atomic
      d2pthinf(iseqp,1) = d2pthinf(iseqp,1) - d1pthflw_out(ipth,ilev)             !! total inflow [m3/s] (from upstream)
    endif
    d1pthflw_pre(ipth,ilev)=d1pthflw_out(ipth,ilev)
  end do

end do
!$omp end parallel do

end subroutine calc_pthout
end module calc_pthout_mod
