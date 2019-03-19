SUBROUTINE CaMaMAIN(iyyyy,imm,idd,istep,r2roffin)

use parkind1  ,only: jpim, jprm
use mod_input ,only: dt, nxin, nyin
use mod_time  ,only: kstep, kminstart, &
                     kyyyymmdd,  khhmm,  kmin, &
                     kyyyymmddp, khhmmp, kminp
use lib_dates ,only: min2date, splitdate

! --------------------------------------------

use control_inp_mod
use control_phy_mod
use control_out_mod
use control_rest_mod
use omp_lib

! --------------------------------------------

IMPLICIT NONE

integer(kind=jpim), intent(in) :: istep               !! local time step
real(kind=jprm), intent(in) :: r2roffin(nxin,nyin)    !! input runoff (mm/day)
integer(kind=jpim), intent(inout) :: iyyyy, imm, idd  !! date
integer(kind=jpim) :: iyyyyp,immp,iddp                !! date (previous)

print*,'[1] CaMa Main '
! ==========================================
! 1. update time 
  iyyyyp=iyyyy
  immp=imm
  iddp=idd
  kyyyymmddp=kyyyymmdd
  khhmmp=khhmm
  kminp=kmin
  kstep=istep
  kmin=int ( kminstart+kstep*(dt/60_jpim) )
  call min2date(kmin,kyyyymmdd,khhmm)
  call splitdate(kyyyymmdd,iyyyy,imm,idd)

! ==========================================
! 2. update runoff at the CoLM timestep
  call control_inp(r2roffin)

! ==========================================
! 3. main calculation
  call control_phy
  
! ============================================
! 4. update diagnostics, and write to output if thats the case
  call control_out(iyyyy,immp,iddp)
                   
!  ===============================
! 5. wrte_restart
  call control_rest(iyyyy,imm,idd)

!!==================================================

END SUBROUTINE CaMaMAIN
