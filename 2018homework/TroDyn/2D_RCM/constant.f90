!define constants
!by renql 20190412

module constant
implicit none

public
    integer :: i, ilat !used to assignment or cycling
    real    :: f
    integer, parameter :: pr = 8
    integer, parameter :: nday = 365   
    integer, parameter :: nlat = 23     !latitude spacing is 8, from -88 to 88
    integer, parameter :: nlev = 10     !latitude spacing is 8, from -88 to 88
    real(kind=pr), parameter, dimension(nday) :: day = (/(i,i=1,nday)/)
    real(kind=pr), parameter, dimension(nlat) :: lat = (/(f,f=-88,88,8)/)

    real(kind=pr), parameter :: cp = 1004 !J/(kg*K),specific heat at constant pressure for air
    
contains

    subroutine calc_density(p,t,density)
        implicit none
        real(kind=pr), intent(in) :: p,t !presure(Pa),temperature(K)
        real(kind=pr), intent(out):: density !kg/m3, density for air
        real(kind=pr), parameter :: p0 = 100000 !Pa, standard atmospheric pressure
        real(kind=pr), parameter :: t0 = 273.15 !K, standard temperature

        density = 1.293*(p/p0)*(t0*t)
    end subroutine calc_density

end module constant 
