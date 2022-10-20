SUBROUTINE CaMaINI(iyyyy,imm,idd,eyyyy,emm,edd)

use parkind1  ,only: jpim
use mod_input ,only: lognam, clogout, isyear, &
                     isyear, ismon, isday, ieyear, iemon, ieday
use mod_map   ,only: regionall, regionthis

use mod_time  ,only: kyyyymmdd, khhmm, kmin, &
                     kyyyymmddp, khhmmp, kminp, nsteps
use lib_dates ,only: splitdate

! --------------------------------------------

use init_inputnam_mod
use init_map_mod
use init_topo_mod
use init_cond_mod

use omp_lib

! --------------------------------------------

IMPLICIT NONE

integer(kind=jpim), intent(inout) :: iyyyy, imm, idd   !! date
integer(kind=jpim), intent(in)    :: eyyyy, emm, edd   !! date

print*, '[1] CaMa INI'
!!================================================
! *** 1. read input namelist 
call init_inputnam

isyear = iyyyy; ismon  = imm; isday  = idd
ieyear = eyyyy; iemon  = emm; ieday  = edd

if (regionall>=2 )then   !! for mpi
  write(clogout,'(a5,i4.4,a1,i2.2,a4)') './out', isyear, '-', regionthis, '.log'
  lognam=12
  open(lognam,file=clogout,form='formatted')
endif

! *** 2. read input river map
call init_map

! *** 3. set topography
call init_topo

! *** 4. set initial conditions (zero-storage / restart)
call init_cond(iyyyy,imm,idd)

call init_time
call init_diag
call splitdate(kyyyymmdd,iyyyy,imm,idd)

!!==================================================

END SUBROUTINE CaMaINI
! ------ EOP ---------
