module calc_fldstg_mod
contains
subroutine calc_fldstg
! ================================================
! calculate river and floodplain staging
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: nlfp
use mod_map    ,only: nseqall
use mod_map    ,only: d2grarea, d2rivlen, d2rivwth, d2rivstomax, d2fldstomax, d2fldgrd, dfrcinc
use mod_prog   ,only: d2rivsto_out, d2fldsto_out
use mod_diag   ,only: d2rivdph, d2flddph, d2fldfrc, d2fldare, dglbrivsto, dglbfldsto, dglbfldare
implicit none

!*** local
!$ save
integer(kind=jpim) ::  iseq, i, ii
real(kind=jprb)    ::  dstoall, dstonow, dstopre, dwthnow, dwthpre, ddphpre, dwthinc
!$omp threadprivate (i,dstoall, dstonow, dstopre, dwthnow, dwthpre, ddphpre, dwthinc)
!!==============================
dglbrivsto=0.d0
dglbfldsto=0.d0
dglbfldare=0.d0

!$omp parallel do reduction(+:dglbrivsto,dglbfldsto,dglbfldare)
do iseq=1, nseqall
!
  dstoall = d2rivsto_out(iseq,1) + d2fldsto_out(iseq,1)
  if( dstoall > d2rivstomax(iseq,1) )then
    i=1
    dstopre = d2rivstomax(iseq,1)
    dwthpre = d2rivwth(iseq,1)
    ddphpre = 0.d0
    dwthinc = d2grarea(iseq,1) * d2rivlen(iseq,1)**(-1.) * dfrcinc
!*  do while( dstoall > d2fldstomax(iseq,1,i) .and. i<=nlfp)
    do ii = 1, nlfp                              ! /04/2014/ yongjiu dai
       if( dstoall > d2fldstomax(iseq,1,ii) )then
         dstopre = d2fldstomax(iseq,1,ii)
         dwthpre = dwthpre + dwthinc
         ddphpre = ddphpre + d2fldgrd(iseq,1,ii) * dwthinc
         i=i+1
!*       if( i>nlfp ) exit
       end if
    end do
    if( i>nlfp )then
      dstonow = dstoall - dstopre
      dwthnow = 0.d0
      d2flddph(iseq,1) = ddphpre + dstonow * dwthpre**(-1.) * d2rivlen(iseq,1)**(-1.)
    else
      dstonow =  dstoall - dstopre
      dwthnow = -dwthpre + ( dwthpre**2. + 2.d0 * dstonow * d2rivlen(iseq,1)**(-1.) * d2fldgrd(iseq,1,i)**(-1.) )**0.5
      d2flddph(iseq,1) = ddphpre + d2fldgrd(iseq,1,i) * dwthnow
    endif
    d2rivsto_out(iseq,1) = d2rivstomax(iseq,1) + d2rivlen(iseq,1) * d2rivwth(iseq,1) * d2flddph(iseq,1)
    d2rivdph(iseq,1) = d2rivsto_out(iseq,1) * d2rivlen(iseq,1)**(-1.) * d2rivwth(iseq,1)**(-1.)
!
    d2fldsto_out(iseq,1) = dstoall - d2rivsto_out(iseq,1)
    d2fldsto_out(iseq,1) = max( d2fldsto_out(iseq,1), 0.d0 )
    d2fldfrc(iseq,1) = (-d2rivwth(iseq,1) + dwthpre + dwthnow ) * (dwthinc*10.d0)**(-1.)
    d2fldfrc(iseq,1) = max( d2fldfrc(iseq,1),0.d0)
    d2fldfrc(iseq,1) = min( d2fldfrc(iseq,1),1.d0)
    d2fldare(iseq,1) = d2grarea(iseq,1)*d2fldfrc(iseq,1)
  else
    d2rivsto_out(iseq,1) = dstoall
    d2rivdph(iseq,1) = dstoall * d2rivlen(iseq,1)**(-1.) * d2rivwth(iseq,1)**(-1.)
    d2rivdph(iseq,1) = max( d2rivdph(iseq,1), 0.d0 )
    d2fldsto_out(iseq,1) = 0.d0
    d2flddph(iseq,1) = 0.d0
    d2fldfrc(iseq,1) = 0.d0
    d2fldare(iseq,1) = 0.d0
  endif

  dglbrivsto      = dglbrivsto + d2rivsto_out(iseq,1)
  dglbfldsto      = dglbfldsto + d2fldsto_out(iseq,1)
  dglbfldare      = dglbfldare + d2fldare(iseq,1)

end do
!$omp end parallel do


end subroutine calc_fldstg
end module calc_fldstg_mod
