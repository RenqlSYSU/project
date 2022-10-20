module calc_damout_mod
contains
subroutine calc_damout
! ================================================
! calculate discharge, inertial equation
! ================================================
use parkind1   ,only: jpim, jprb
use mod_map    ,only: i1next,   i2vector
use mod_prog   ,only: d2rivout_out, d2fldout_out
use mod_diag   ,only: d2rivinf,     d2fldinf

implicit none

!*** local
      integer(kind=jpim) ::  ix, iy, iseq, jseq
      real(kind=jprb)    ::  doutflw, dexcess, drivred, dfldred       !! total outflw, excess flow, river flow reduction, floodplain flow reduction

!*** dam parameter
      real(kind=jprb)    ::  dmaxout
! ================================================

!! for stung toreng (e105,94,n13.55)

ix=160   !! please check (ix,iy) of the target point
iy=215
iseq=i2vector(ix,iy)  ! convert 2d-map (ix,iy) to 1d-vector (iseq,1) 
jseq=i1next(iseq)

!! dam operation (modify rivout & fldout)

doutflw=d2rivout_out(iseq,1)+d2fldout_out(iseq,1)
dmaxout=40000.

if( doutflw > dmaxout )then
  dexcess=doutflw-dmaxout
  if( d2fldout_out(iseq,1)>dexcess )then
    d2fldout_out(iseq,1)=d2fldout_out(1,iseq)-dexcess
    d2fldinf(jseq,1)=d2fldinf(jseq,1)-dexcess
  else
    drivred=d2rivout_out(iseq,1)-dmaxout
    dfldred=d2fldout_out(iseq,1)

    d2rivout_out(iseq,1)=dmaxout
    d2fldout_out(iseq,1)=0.0

    d2rivinf(jseq,1)=d2rivinf(jseq,1)-drivred
    d2fldinf(jseq,1)=d2fldinf(jseq,1)-dfldred
  endif
endif


end subroutine calc_damout
end module calc_damout_mod
