!calculate irradiance at spectral bands (W/(m2*um)) according Blackbody radiation theorem
!calculate the sum of irradiance along all spectral
!by renql 20190413

module spectral
implicit none

contains
    subroutine swbnd_irradiance(temperature)
        use constant, only: nbndsw,swbnd

        implicit none
        real(kind=pr), intent(in) :: p,t !presure(Pa),temperature(K)
        real(kind=pr), intent(out):: density !kg/m3, density for air
        real(kind=pr), parameter :: p0 = 100000 !Pa, standard atmospheric pressure
        real(kind=pr), parameter :: t0 = 273.15 !K, standard temperature

        density = 1.293*(p/p0)*(t0*t)
    end subroutine calc_density

