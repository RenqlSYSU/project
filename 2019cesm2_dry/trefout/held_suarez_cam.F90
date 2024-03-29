#define MODHS 1
#undef MODHS
! ---------------------------
! Modified to output an additional history field TREF which contains
! the relaxation temperature profile
! Isla Simpson 3rd June 2017
! 
! Modifications are denoted by 
!
! !OUTFLD
! blah blah blah
! !END-OUTFLD
!
! Updated for CESM2 release 30th May 2018 (Isla Simpson)
! ----------------------------
module held_suarez_cam
  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pcols, pver

  implicit none
  private
  save

  public :: held_suarez_init, held_suarez_tend

  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

!======================================================================= 
contains
!======================================================================= 

  subroutine held_suarez_init(pbuf2d)
    use physics_buffer,     only: physics_buffer_desc
    use cam_history,        only: addfld, add_default
    use physconst,          only: cappa, cpair
    use ref_pres,           only: pref_mid_norm, psurf_ref
    use held_suarez,        only: held_suarez_1994_init

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    ! Set model constant values
    call held_suarez_1994_init(cappa, cpair, psurf_ref, pref_mid_norm)

    ! This field is added by radiation when full physics is used
    call addfld('QRS', (/ 'lev' /), 'A', 'K/s', &
         'Temperature tendency associated with the relaxation toward the equilibrium temperature profile')

    !OUTFLD
    call addfld('TREF', (/ 'lev' /), 'A', 'K', &
         'Equilibrium temperature profile')
    !END-OUTFLD

    call add_default('QRS', 1, ' ')
    !OUTFLD
    call add_default('TREF', 1, ' ')
    !END-OUTFLD

 end subroutine held_suarez_init

  subroutine held_suarez_tend(state, ptend, ztodt)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !  algorithm 1: Held/Suarez IDEALIZED physics
    !  algorithm 2: Held/Suarez IDEALIZED physics (Williamson modified stratosphere
    !  algorithm 3: Held/Suarez IDEALIZED physics (Lin/Williamson modified strato/meso-sphere
    !
    ! Author: J. Olson
    ! 
    !-----------------------------------------------------------------------
    use physconst,          only: cpairv
    use phys_grid,          only: get_rlat_all_p
    use physics_types,      only: physics_state, physics_ptend
    use physics_types,      only: physics_ptend_init
    use cam_abortutils,     only: endrun
    use cam_history,        only: outfld
    use held_suarez,        only: held_suarez_1994

    !
    ! Input arguments
    !
    type(physics_state), intent(inout) :: state
    real(r8),            intent(in)    :: ztodt            ! Two times model timestep (2 delta-t)
                                                           !
                                                           ! Output argument
                                                           !
    type(physics_ptend), intent(out)   :: ptend            ! Package tendencies
                                                           !
    !---------------------------Local workspace-----------------------------

    integer                            :: lchnk            ! chunk identifier
    integer                            :: ncol             ! number of atmospheric columns

    real(r8)                           :: clat(pcols)      ! latitudes(radians) for columns
    real(r8)                           :: pmid(pcols,pver) ! mid-point pressure
    integer                            :: i, k             ! Longitude, level indices

    !OUTFLD
    real(r8)                           :: trefout(pcols,pver) !Equilibrium temperature profile for output
    !END-OUTFLD

    !
    !-----------------------------------------------------------------------
    !

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_rlat_all_p(lchnk, ncol, clat)
    do k = 1, pver
      do i = 1, ncol
        pmid(i,k) = state%pmid(i,k)
      end do
    end do

    ! initialize individual parameterization tendencies
    call physics_ptend_init(ptend, state%psetcols, 'held_suarez', ls=.true., lu=.true., lv=.true.)

    ! OUTFLD
    ! Modify the call to held_suarez_1994 to return the relaxation temperature
    ! profile (trefout)
    !call held_suarez_1994(pcols, ncol, clat, state%pmid, &
    !     state%u, state%v, state%t, ptend%u, ptend%v, ptend%s)
     call held_suarez_1994(pcols, ncol, clat, pmid, &
         state%u, state%v, state%t, ptend%u, ptend%v, ptend%s, trefout)
    ! END-OUTFLD

    ! Note, we assume that there are no subcolumns in simple physics
    pmid(:ncol,:) = ptend%s(:ncol, :)/cpairv(:ncol,:,lchnk)
    if (pcols > ncol) then
      pmid(ncol+1:,:) = 0.0_r8
    end if
    call outfld('QRS', pmid, pcols, lchnk)
    !OUTFLD
    call outfld('TREF', trefout, pcols, lchnk)
    !END-OUTFLD

  end subroutine held_suarez_tend

end module held_suarez_cam
