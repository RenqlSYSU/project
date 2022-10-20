module calc_outpre_mod
contains
subroutine calc_outpre
! ================================================
! to   calculate discharge, diffusive wave, for storage only restart
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: lfldout, lpthout
use mod_input  ,only: pmanriv,  pmanfld,  pdstmth
use mod_map    ,only: i1next,   nseqall,  nseqriv,  nseqmax,  npthout
use mod_map    ,only: d2rivelv, d2elevtn, d2nxtdst, d2rivwth, d2rivlen
use mod_map    ,only: npthout, npthlev, pth_upst, pth_down, pth_dst, pth_elv, pth_wth, pth_man
use mod_prog   ,only: d2rivout_pre, d2fldout_pre, d1pthflw_pre ,d2rivdph_pre         !! output
use mod_prog   ,only: d2fldsto_out                                    !! input
use mod_diag   ,only: d2rivdph, d2sfcelv, d2flddph                    !! input

implicit none

!*** local
      real(kind=jprb)    ::  d2sfcelv_fld(nseqmax,1)                  !! water surface elev (t)  [m] (for floodplain)
      real(kind=jprb)    ::  d2flware_now(nseqmax,1)                  !! flow area (t)           [m2]

!$ save
      integer(kind=jpim) ::  iseq, jseq, ipth, ilev, iseqp, jseqp
      real(kind=jprb)    ::  dslope, darea, dflw
!$omp threadprivate         (dslope, darea, dflw, jseq, ilev, iseqp, jseqp)
! ================================================
!$omp parallel do
do iseq=1, nseqall
  d2rivdph_pre(iseq,1) = d2rivdph(iseq,1)                           !! bugfix v362
  d2sfcelv(iseq,1)     = d2rivelv(iseq,1) + d2rivdph(iseq,1)
  
  d2sfcelv_fld(iseq,1) = max( d2elevtn(iseq,1), d2sfcelv(iseq,1) )
  d2flware_now(iseq,1) = d2fldsto_out(iseq,1) * d2rivlen(iseq,1)**(-1.d0)
end do
!$omp end parallel do

!$omp parallel do
do iseq=1, nseqriv                                                    !! for normal cells
  jseq=i1next(iseq)

  dslope =  (d2sfcelv(iseq,1)-d2sfcelv(jseq,1)) * d2nxtdst(iseq,1)**(-1.d0)
  dflw = max(d2sfcelv(iseq,1),d2sfcelv(jseq,1)) - d2rivelv(iseq,1)
  darea = d2rivwth(iseq,1) * dflw

  if( darea>1.d-5 )then
    d2rivout_pre(iseq,1) = darea * ( pmanriv**(-1.d0) * dflw**0.6667d0 * abs(dslope)**(0.5d0) )
    if( dslope<0.d0 ) d2rivout_pre(iseq,1)=-d2rivout_pre(iseq,1)
  else
    d2rivout_pre(iseq,1) = 0.d0
  endif
end do
!$omp end parallel do

!$omp parallel do                                                     !! for river mouth grids
do iseq=nseqriv+1, nseqall
  dslope = ( d2sfcelv(iseq,1)-d2elevtn(iseq,1) ) * pdstmth**(-1.d0)
!!  dslope = max( dslope, 0.d0 )                                        !! no backfrom from ocean
  dflw   = d2rivdph(iseq,1)
  darea  = d2rivwth(iseq,1) * dflw
  if( darea>1.d-5 )then
    d2rivout_pre(iseq,1) = darea * ( pmanriv**(-1.d0) * dflw**0.6667d0 * abs(dslope)**(0.5d0) )
    if( dslope<0.d0 ) d2rivout_pre(iseq,1)=-d2rivout_pre(iseq,1)
  else
    d2rivout_pre(iseq,1) = 0.d0
  endif
end do
!$omp end parallel do



if( lfldout )then
 !$omp parallel do
 do iseq=1, nseqriv                                                    !! for normal cell
  jseq=i1next(iseq)

  dslope  = (d2sfcelv_fld(iseq,1)-d2sfcelv_fld(jseq,1)) * d2nxtdst(iseq,1)**(-1.d0)

  dflw = min(d2sfcelv_fld(iseq,1),d2sfcelv_fld(jseq,1)) - d2elevtn(iseq,1)        !! should be min[i,j] to avoid instability
  dflw = max(dflw,0.d0)
  darea= min( d2flware_now(iseq,1), d2flware_now(jseq,1) )
  if( darea<=0.d0 ) darea=d2flware_now(iseq,1)

  if( darea>1.d-5 )then
    d2fldout_pre(iseq,1) = darea * ( pmanfld**(-1.d0) * dflw**0.6667d0 * abs(dslope)**(0.5d0) )
    if( dslope<0.d0 ) d2fldout_pre(iseq,1)=-d2fldout_pre(iseq,1)
  else
    d2fldout_pre(iseq,1) = 0.d0
  endif
end do
!$omp end parallel do

!$omp parallel do
do iseq=nseqriv+1, nseqall                                            !! for river mouth cell
  dslope = ( d2sfcelv_fld(iseq,1)-d2elevtn(iseq,1) ) * pdstmth**(-1.d0)
!  dslope = max( dslope, 0.d0 )                                        !! no backward flow from ocean
  dflw   = d2flddph(iseq,1)

  darea= d2flware_now(iseq,1)

  if( darea>1.d-5 )then
    d2fldout_pre(iseq,1) = darea * ( pmanfld**(-1.d0) * dflw**0.6667d0 * abs(dslope)**(0.5d0) )
    if( dslope<0.d0 ) d2fldout_pre(iseq,1)=-d2fldout_pre(iseq,1)
  else
    d2fldout_pre(iseq,1) = 0.d0
  endif
 end do
 !$omp end parallel do
endif



if( lpthout )then
 !$omp parallel do
 do ipth=1, npthout  
  iseqp=pth_upst(ipth)
  jseqp=pth_down(ipth)

  dslope  = (d2sfcelv(iseqp,1)-d2sfcelv(jseqp,1)) * pth_dst(ipth)**(-1.d0)

  do ilev=1, npthlev
    dflw = max(d2sfcelv(iseqp,1),d2sfcelv(jseqp,1)) - pth_elv(ipth,ilev) 
    dflw = max(dflw,0.d0)

    if( dflw>1.d-5 )then
      d1pthflw_pre(ipth,ilev) = pth_wth(ipth,ilev) * dflw * ( pth_man(ilev)**(-1.d0) * dflw**0.6667d0 * abs(dslope)**(0.5d0) )
      if( dslope<0.d0 ) d1pthflw_pre(ipth,ilev)=-d1pthflw_pre(ipth,ilev)
    else
      d1pthflw_pre(ipth,ilev) = 0.d0
    endif
  end do
 end do
 !$omp end parallel do
endif

end subroutine calc_outpre
end module calc_outpre_mod
