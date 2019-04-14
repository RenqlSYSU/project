!define constants
!by renql 20190412

module constant
implicit none

private
    integer :: i, ilat !used to assignment or cycling
    real    :: f

public
    integer, parameter :: pr = 8
    integer, parameter :: nday = 365   
    integer, parameter :: nlat = 23     !latitude spacing is 8, from -88 to 88
    integer, parameter :: nlev = 10     !latitude spacing is 8, from -88 to 88
    integer, parameter :: nbndsw = 17   !number of shorwave spectral intervals,from 0.2 to 5 micro-meters
    integer, parameter :: nbndlw = 17   !number of shorwave spectral intervals,from 0.2 to 5 micro-meters

    real(kind=pr), parameter, dimension(nday) :: day = (/(i,i=1,nday)/)
    real(kind=pr), parameter, dimension(nlat) :: lat = (/(f,f=-88,88,8)/)
    real(kind=pr), parameter, dimension(nbndsw) :: sw_lgth = (/(f,f=0.2,5,0.3)/)
    real(kind=pr), parameter, dimension(nbndlw) :: lw_lgth = (/(f,f=0.2,5,0.3)/)
    real(kind=pr), parameter, dimension(nbndsw) :: swbnd = 0.3
    real(kind=pr), parameter, dimension(nbndlw) :: lwbnd = (/(f,f=0.2,5,0.3)/)

    real(kind=pr), parameter :: cp = 1004 !J/(kg*K),specific heat at constant pressure for air
    real(kind=pr), parameter :: pai = 3.1415926 
    real(kind=pr), parameter :: Plank_cons = 6.626*(10**(-34)) !Plank constant (J/s) 
    real(kind=pr), parameter :: Boltz_cont = 1.38065*(10**(-23)) !Boltzmann constant, J/k
    real(kind=pr), parameter :: vlight = 3*(10**8) !velocity of light, m/s
    
contains

    subroutine calc_density(p,temp,density)
        implicit none
        real(kind=pr), intent(in) :: p,temp !presure(Pa),temperature(K)
        real(kind=pr), intent(out):: density !kg/m3, density for air
        real(kind=pr), parameter :: p0 = 100000 !Pa, standard atmospheric pressure
        real(kind=pr), parameter :: t0 = 273.15 !K, standard temperature

        density = 1.293*(p/p0)*(t0*t)
    end subroutine calc_density

    subroutine calc_bnd_irradiance(temp,nbnd,lgth,irrad)
        integer, intent(in) :: nbnd !temperature(K)
        real(kind=pr), intent(in) :: temp !temperature(K)
        real(kind=pr), dimension(nbnd),intent(in):: lgth !um
        real(kind=pr), dimension(nbnd),intent(out):: irrad !W/(m2*um)

        irrad = 2*pai*Plank_cons*(vlight**2)/(lgth**5)&
                &/(exp(vlight*Plank_cons/(Boltz_cont*lgth*temp))-1)
    end subroutine calc_bnd_irradiance

    subroutine calc_all_irradiance(irrad,nbnd,bnd,radi)
        integer, intent(in) :: nbnd !temperature(K)
        real(kind=pr), dimension(nbnd),intent(in):: irrad !W/(m2*um)
        real(kind=pr), dimension(nbnd),intent(in):: bnd !W/(m2*um)
        real(kind=pr), intent(out):: radi !W/m2
        
        radi = sum(irrad*bnd)
    end subroutine calc_all_irradiance

end module constant 
