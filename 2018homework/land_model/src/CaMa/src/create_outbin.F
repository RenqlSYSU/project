module create_outbin_mod

contains

subroutine create_outbin(iyear)

!**** *create_outbin*  - create_out binary files

!     --------------
!        original : 2011-05-04
!     ------------------------------------------------------------------

use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam,loutvec,nx,ny,lpthout
use mod_map    ,only: nseqmax, npthout, npthlev
use mod_output ,only: nvars,funits,coutvardir,coutvn,lopen,&
                      csufbin, csufvec, csufpth
implicit none
integer             ::  iyear               !! year
integer             ::  jf

!!=== 
if (lopen) then
  write(lognam,*) "create_outbin: closing previous units"

  do jf=1, nvars
    if (coutvardir(jf).ne.'NONE') close(funits(jf))
  enddo
endif

!!====
write(lognam,*) "****create_outbin: ****"

if (loutvec) then ! vectorial output
  write(lognam,*) "create_outbin: vectorial output"

  do jf=1, nvars
    if (coutvardir(jf).ne.'NONE')then
      if( jf/=nvars) then
        call open_outvec(iyear,funits(jf), coutvardir(jf),coutvn(jf))
      else
        call open_outpth(iyear,funits(jf), coutvardir(jf),coutvn(jf))
      endif
    endif
  end do
  lopen=.true.

else ! nx,ny output
  write(lognam,*) "create_outbin: regular nx,ny ouput"
  do jf=1, nvars
    if (coutvardir(jf).ne.'NONE')then
      if( jf/=nvars) then
        call open_outbin(iyear,funits(jf), coutvardir(jf),coutvn(jf))
      else
        call open_outpth(iyear,funits(jf), coutvardir(jf),coutvn(jf))
      endif
    endif
  end do
  
  lopen=.true.
endif

contains


!!================================================
subroutine open_outvec(iyear,ifn,coutdir,cvarname)
implicit none

! 
integer             ::  iyear               !! year
integer             ::  ifn                 !! file num
character*128       ::  coutdir             !! output dir
character*128       ::  cvarname            !! variavle name
! 
character*128       ::  cfile
! 
integer             ::  len_trim
integer             ::  ilen1, ilen2
character*128       ::  clen1, clen2
!!==============================
ilen1=len_trim(coutdir)
write(clen1,*) ilen1
ilen2=len_trim(cvarname)
write(clen2,*) ilen2
!
write(cfile,'(a'//clen1//',a'//clen2//',i4,a10)') coutdir,cvarname,iyear,csufvec
open(ifn,file=cfile,form='unformatted',access='direct',recl=4*nseqmax)

end subroutine open_outvec





!!================================================
subroutine open_outpth(iyear,ifn,coutdir,cvarname)
implicit none

! 
integer             ::  iyear               !! year
integer             ::  ifn                 !! file num
character*128       ::  coutdir             !! output dir
character*128       ::  cvarname            !! variavle name
! 
character*128       ::  cfile
! 
integer             ::  len_trim
integer             ::  ilen1, ilen2
character*128       ::  clen1, clen2
!!==============================
ilen1=len_trim(coutdir)
write(clen1,*) ilen1
ilen2=len_trim(cvarname)
write(clen2,*) ilen2
!
write(cfile,'(a'//clen1//',a'//clen2//',i4,a10)') coutdir,cvarname,iyear,csufpth
open(ifn,file=cfile,form='unformatted',access='direct',recl=4*npthout*npthlev)

end subroutine open_outpth





!!==================================================
subroutine open_outbin(iyear,ifn,coutdir,cvarname)
implicit none

integer             ::  iyear               !! year
integer             ::  ifn                 !! file num
character*128       ::  coutdir             !! output dir
character*128       ::  cvarname            !! variavle name

character*128       ::  cfile

integer             ::  len_trim
integer             ::  ilen1, ilen2
character*128       ::  clen1, clen2
!!==============================
ilen1=len_trim(coutdir)
write(clen1,*) ilen1
ilen2=len_trim(cvarname)
write(clen2,*) ilen2

write(cfile,'(a'//clen1//',a'//clen2//',i4,a10)') coutdir,cvarname,iyear,csufbin
open(ifn,file=cfile,form='unformatted',access='direct',recl=4*nx*ny)

end subroutine open_outbin

end subroutine create_outbin
end module create_outbin_mod
