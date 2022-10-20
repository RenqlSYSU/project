module calc_watbal_mod
contains      
subroutine calc_watbal(it)
! ================================================
! write water balance monitoring to iofile
! ================================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: lognam,dt,dt_def
use mod_time   ,only: kminstart,kstep
use mod_diag   ,only: dglbrivsto,dglbfldsto,dglbstopre,dglbstonxt,dglbstonew,dglbrivinf,dglbrivout,dglbfldare
use lib_dates  ,only: splitdate,min2date

implicit none

integer(kind=jpim),intent(in) :: it        !! step in adaptive time loop

!*** local
real(kind=jprb)    ::  dmissing            !! water ballance error [kg]
integer(kind=jpim) ::  iyear, imon, iday, ihour, imin,jkmin,jkyyyymmdd,jkhhmm

!*** parameter
real(kind=jprb)    ::  dord
parameter             (dord=1.d-9)
! ================================================
! write water balance
! ================================================

jkmin=int ( kminstart+(kstep-1)*(dt_def/60_jpim) + it*dt/60_jprb )
call min2date(jkmin,jkyyyymmdd,jkhhmm)
call splitdate(jkyyyymmdd,iyear,imon,iday)
ihour=jkhhmm/100
imin=jkhhmm-ihour*100


dmissing   = dglbstopre - dglbstonxt + dglbrivinf - dglbrivout
write(lognam,'(i4.4,4(a1,i2.2),i6,3f12.3,2x,2f12.3,2x,2f12.3,g12.3,f12.3)') &
  iyear, '/', imon, '/', iday, '_', ihour, ':', imin, it,   &
  dglbstopre*dord, dglbstonxt*dord, dglbstonew*dord ,dglbrivsto*dord, dglbfldsto*dord, &
  dglbrivinf*dord, dglbrivout*dord, dmissing*dord,   dglbfldare*dord

end subroutine calc_watbal
end module calc_watbal_mod
