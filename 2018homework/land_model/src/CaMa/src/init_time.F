subroutine init_time
!**** *init_cond*  - init time

!     --------------
!        original : 2011-05-04
!     ------------------------------------------------------------------

use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam,&
                      isyear,ismon,isday,&
                      ieyear,iemon,ieday,&
                      dt
use mod_time   ,only: kstep,nsteps,kyyyymmdd,khhmm,kmin,kminstart

use lib_dates  ,only: date2min
implicit none

integer(kind=jpim)     :: syyyymmdd,shhmm,eyyyymmdd,ehhmm
integer(kind=jpim)     :: emin

!!===========================
write(lognam,*) " "
write(lognam,*) "****init_time: ****"

kstep=0 

syyyymmdd=isyear*10000+ismon*100+isday
shhmm=0_jpim
eyyyymmdd=ieyear*10000+iemon*100+ieday
ehhmm=0_jpim

kmin=date2min(syyyymmdd,shhmm)
emin=date2min(eyyyymmdd,ehhmm)
nsteps=int ( ( (emin-kmin)*60_jpim )/dt )

write(lognam,*) 'start date:',syyyymmdd,shhmm,kmin
write(lognam,*) 'end date:',eyyyymmdd,ehhmm,emin
write(lognam,*) 'nsteps:',nsteps

kyyyymmdd=syyyymmdd
khhmm=shhmm
kminstart=kmin


end subroutine init_time
