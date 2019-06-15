module held_suarez
! ----------------------------------------------------------------------------------
! Modified to inlude an option for the Polvani and Kushner (2002) , 
! GRL, 29, 7, 10.1029/2001GL014284 (PK02) equilibrium temperature profile.
! 
! If pkstrat = .True. then Polvani and Kushner relaxation temperature profile
! is used.  
!
! Namelist parameter: vgamma, sets the vortex strength (gamma parameter in PK02)
!
! Modifications are denoted by 
! 
! !PKSTRAT
! blah blah blah
! !END-PKSTRAT
!
! Isla Simpson 8th June 2017
!
! Modified for CESM2 release, Isla Simpson, 30th May 2018
! ----------------------------------------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  !PKSTRAT
  use cam_logfile, only: iulog
  !END-PKSRAT


  implicit none
  private
  save

  public :: held_suarez_1994_init
  public :: held_suarez_1994

  !!
  !! Forcing parameters
  !!
  real(r8), parameter :: efoldf  =  1._r8  ! efolding time for wind dissipation
  real(r8), parameter :: efolda  = 40._r8  ! efolding time for T dissipation
  real(r8), parameter :: efolds  =  4._r8  ! efolding time for T dissipation
  real(r8), parameter :: sigmab  =  0.7_r8 ! threshold sigma level
  real(r8), parameter :: t00     = 200._r8 ! minimum reference temperature
  real(r8), parameter :: kf      = 1._r8/(86400._r8*efoldf) ! 1./efolding_time for wind dissipation

  real(r8), parameter :: onemsig = 1._r8 - sigmab ! 1. - sigma_reference

  real(r8), parameter :: ka      = 1._r8/(86400._r8 * efolda) ! 1./efolding_time for temperature diss.
  real(r8), parameter :: ks      = 1._r8/(86400._r8 * efolds)

  !!
  !! Model constants, reset in init call
  !!
  real(r8)              :: cappa = 2.0_r8 / 7.0_r8  ! R/Cp
  real(r8)              :: cpair = 1004.0_r8        ! specific heat of dry air (J/K/kg)
  real(r8)              :: psurf_ref = 0.0_r8       ! Surface pressure
  ! pref_mid_norm are layer midpoints normalized by surface pressure ('eta' coordinate)
  real(r8), allocatable :: pref_mid_norm(:)
  integer               :: pver                     ! Num vertical levels

  !PKSTRAT
  ! Forcing parameters for PK02 option
  integer, parameter :: lat0=-50.
  integer, parameter :: dellat=10.
  integer, parameter :: dely=60.
  integer, parameter :: eps=-10. !(epsilon parameter in A4 of PK)
  integer, parameter :: delz=10.
  real(r8), parameter :: efoldstrat = 0.5_r8
  real(r8), parameter :: kfstrat     = 1._r8/(86400._r8*efoldstrat)
  integer, parameter :: pret=0.5*100. ! lower limit (Pa) for upper level damping
  !END-PKSTRAT


!======================================================================= 
contains
!======================================================================= 

  subroutine held_suarez_1994_init(cappa_in, cpair_in, psurf_ref_in, pref_mid_norm_in)
    !! Dummy arguments
    real(r8), intent(in) :: cappa_in
    real(r8), intent(in) :: cpair_in
    real(r8), intent(in) :: psurf_ref_in
    real(r8), intent(in) :: pref_mid_norm_in(:)

    pver = size(pref_mid_norm_in)
    allocate(pref_mid_norm(pver))
    cappa         = cappa_in
    cpair         = cpair_in
    psurf_ref     = psurf_ref_in
    pref_mid_norm = pref_mid_norm_in

  end subroutine held_suarez_1994_init

!  subroutine held_suarez_1994(pcols, ncol, clat, pmid, &
!       u, v, t, du, dv, s)
!PKSTRAT
   subroutine held_suarez_1994(pcols, ncol, clat, pmid, &
       u, v, t, du, dv, s, pkstrat, vgamma)
!END-PKSTRAT

    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Implement idealized Held-Suarez forcings
    !    Held, I. M., and M. J. Suarez, 1994: 'A proposal for the
    !    intercomparison of the dynamical cores of atmospheric general
    !    circulation models.'
    !    Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    ! 
    ! 7th March 2017 (Isla Simpson) - modified to include the option of the
    ! Polvani and Kushner (2002), GRL, 29, 7, 10.1029/2001GL014284 stratospheric
    ! relaxation temperature profile
    ! 
    ! Modifications denoted by
    !
    !PKSTRAT
    !blah blah blah
    !END-PKSTRAT
    !
    ! Updated for CESM2 release, Isla Simpson, May 30th
    !-----------------------------------------------------------------------

    !PKSTRAT
    use physconst,      only: rair, gravit ! Gas constant and gravity forpkstrat
    !END-PKSTRAT

    !
    ! Input arguments
    !
    integer,  intent(in)  :: pcols            ! Size of column dimension
    integer,  intent(in)  :: ncol             ! Num active columns
    real(r8), intent(in)  :: clat(pcols)      ! latitudes(radians) for columns
    real(r8), intent(in)  :: pmid(pcols,pver) ! mid-point pressure
    real(r8), intent(in)  :: u(pcols,pver)    ! Zonal wind (m/s)
    real(r8), intent(in)  :: v(pcols,pver)    ! Meridional wind (m/s)
    real(r8), intent(in)  :: t(pcols,pver)    ! Temperature (K)
                                              !
                                              ! Output arguments
                                              !
    real(r8), intent(out) :: du(pcols,pver)   ! Zonal wind tend
    real(r8), intent(out) :: dv(pcols,pver)   ! Meridional wind tend
    real(r8), intent(out) :: s(pcols,pver)    ! Heating rate

    !PKSRAT
    logical, intent(in) :: pkstrat !pkstrat=.True. to use the PK02 TEQ
    real(r8), intent(in) :: vgamma !gamma parameter in PK02 (controling vortex strength) 
    !END-PKSTRAT

    !
    !---------------------------Local workspace-----------------------------
    !
    integer  :: i, k          ! Longitude, level indices

    real(r8) :: kv            ! 1./efolding_time (normalized) for wind
    real(r8) :: kt            ! 1./efolding_time for temperature diss.
    real(r8) :: trefa         ! "radiative equilibrium" T
    real(r8) :: trefc         ! used in calc of "radiative equilibrium" T
    real(r8) :: cossq(ncol)   ! coslat**2
    real(r8) :: cossqsq(ncol) ! coslat**4
    real(r8) :: sinsq(ncol)   ! sinlat**2
    real(r8) :: coslat(ncol)  ! cosine(latitude)

    !PKSTRAT
    real(r8) :: sinlat(ncol) !sin(latitude)
    real(r8) :: deglat(ncol) ! Latitude in degrees for columns
    real(r8) :: pi !pi for computing latitude in degrees from clat
    !END-PKSTRAT


    !
    !-----------------------------------------------------------------------
    !

    !PKSRAT
    real(r8) :: vgammaperm !gamma in K/m (as opposed to K/km) 
    real(r8) :: w, fac1, fac2, delt !weight function (Eq A2 in PK02)

    !US standard atmosphere params
      integer, parameter :: nstd=7 ! # of levels specified in US standard atmosphere
      real(r8) :: pstd_norm(nstd) ! pressures of US standard atmosphere,normalized by PS =101325
      real(r8) :: tstd(nstd) ! Temperutre of standard levels
      real(r8) :: lapse(nstd) ! Lapse rate of standard levels 
      real(r8) :: tvalstd ! Standard atmosphere T at current level
      real(r8) :: tpv ! Polar vortex temperature (Eqn 3 of Kushner and Polvani 2004)
      real(r8) :: pbase,lapsebase,tbase ! US standard parameters at base of layer
      integer :: k2
      real(r8) :: tstd100 ! standard atmosphere temperature at 100hPa


      pstd_norm=(/ 1., 0.223361, 0.0540330, 0.00856668, 0.00109456, 0.000660636, 3.90468e-05/)
      tstd=(/ 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65 /)
      lapse=(/ -0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002 /)

      pi=4.*DATAN(1.D0)
    !END-PKSTRAT




    do i = 1, ncol
      coslat (i) = cos(clat(i))
      sinsq  (i) = sin(clat(i))*sin(clat(i))
      cossq  (i) = coslat(i)*coslat(i)
      cossqsq(i) = cossq (i)*cossq (i)

      !PKSTRAT 
      deglat(i) = (clat(i)/pi)*180.
      sinlat(i) = sin(clat(i))
      !END-PKSTRAT
    end do

    !
    !-----------------------------------------------------------------------
    !
    ! Held/Suarez IDEALIZED physics algorithm:
    !
    !   Held, I. M., and M. J. Suarez, 1994: A proposal for the
    !   intercomparison of the dynamical cores of atmospheric general
    !   circulation models.
    !   Bulletin of the Amer. Meteor. Soc., vol. 75, pp. 1825-1830.
    !
    !-----------------------------------------------------------------------
    !
    ! Compute idealized radiative heating rates (as dry static energy)
    !
    !

    !PK-STRAT
    if (pkstrat) then
       tstd100=216.65_r8*(100._r8/110.906_r8)**(-1._r8*rair*0._r8/9.81_r8) !US standard atmosphere at 100hPa 
       vgammaperm=vgamma/1000._r8 ! convert gamma param into K/m

       !Set up US standard atmosphere
       do k=1,pver
         ! Determine US standard atmosphere params at base of layer
         if (pref_mid_norm(k).lt.pstd_norm(nstd)) then
           pbase=pstd_norm(nstd)
           lapsebase=lapse(nstd)
           tbase=tstd(nstd)
         else
           do k2=1,nstd-1
            if (pref_mid_norm(k).gt.pstd_norm(k2+1)) then
             pbase=pstd_norm(k2)
             lapsebase=lapse(k2)
             tbase=tstd(k2)
             exit
            endif
           end do
         end if

         !US standard atmosphere at level  
         tvalstd=tbase*(pref_mid_norm(k)/pbase)**(-1.*(rair*lapsebase/gravit))

         do i = 1,ncol
           if (pmid(i,k).lt.1e4_r8) then !pre < 100hPa 
             tpv=tstd100*(pmid(i,k)/1e4_r8)**(rair*vgammaperm/gravit)
             w=0.5*(1-tanh( (deglat(i)-lat0)/dellat))
             trefa=(1-w)*tvalstd + w*tpv
           else
             fac1=tstd100
             delt=dely*sinsq(i) + eps*sinlat(i) + delz*log(pmid(i,k)/1e5_r8)*cossq(i)
             fac2=(315._r8-delt)*(pmid(i,k)/1e5_r8)**cappa
             if (fac1.gt.fac2) then
               trefa=fac1
             else
               trefa=fac2
             end if
           end if !end pre < 100hPa

           if (pref_mid_norm(k) > sigmab) then
             kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
             s(i,k) = (trefa - t(i,k))*kt*cpair
           else
             s(i,k) = (trefa - t(i,k))*ka*cpair
           end if

         end do
       end do
     else  !using Held-Suarez TEQ
    !END-PKSTRAT



    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        do i = 1, ncol
          kt = ka + (ks - ka)*cossqsq(i)*(pref_mid_norm(k) - sigmab)/onemsig
          trefc   = 315._r8 - (60._r8 * sinsq(i))
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*kt*cpair
        end do
      else
        do i = 1, ncol
          trefc   = 315._r8 - 60._r8*sinsq(i)
          trefa = (trefc - 10._r8*cossq(i)*log((pmid(i,k)/psurf_ref)))*(pmid(i,k)/psurf_ref)**cappa
          trefa    = max(t00,trefa)
          s(i,k) = (trefa - t(i,k))*ka*cpair
        end do
      end if
    end do

end if !end if pkstrat


    !
    ! Add diffusion near the surface for the wind fields
    !
    do k = 1, pver
      do i = 1, pcols
        du(i,k) = 0._r8
        dv(i,k) = 0._r8
      end do
    end do

    !
    do k = 1, pver
      if (pref_mid_norm(k) > sigmab) then
        kv  = kf*(pref_mid_norm(k) - sigmab)/onemsig
        do i = 1, ncol
          du(i,k) = -kv*u(i,k)
          dv(i,k) = -kv*v(i,k)
        end do
      end if

      !PKSTRAT
      !Add sponge layer in upper levels (see Appendix of PK02
      if (pkstrat) then
        do i = 1,ncol
          if (pmid(i,k).lt.pret) then
            kv=kfstrat*((pret-pmid(i,k))/pret)**(2.)
            du(i,k) = -kv*u(i,k)
            dv(i,k) = -kv*v(i,k)
           endif !pmid < 0.5
        end do
      endif !pkstrat
      !END-PKSTRAT 
    end do

  end subroutine held_suarez_1994

end module held_suarez
