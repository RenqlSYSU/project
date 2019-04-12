!This module 
!by renql 20190412
!

subroutine rad_swf(net_flux)
    implicit none
    
    integer :: i, ilat !used to assignment or cycling
    real    :: f
    integer, parameter :: pr = 8
    integer, parameter :: nday = 365   
    integer, parameter :: nlat = 23     !latitude spacing is 8, from -88 to 88
    integer, parameter :: nlev = 10     !latitude spacing is 8, from -88 to 88
    real(kind=pr), parameter, dimension(nday) :: day = (/(i,i=1,nday)/)
    real(kind=pr), parameter, dimension(nlat) :: lat = (/(f,f=-88,88,8)/)

    real(kind=pr), dimension(nday) :: day_angle
    real(kind=pr), dimension(nday) :: dist_ratio
    real(kind=pr), dimension(nday) :: sun_point
    real(kind=pr), dimension(nday) :: hour_angle
    real(kind=pr), dimension(nday,nlat) :: all_sw !total daily average solar radiation (W/m2)
    real(kind=pr), dimension(nday,nlat) :: zenith !mean zenith angle 

    real(kind=pr), parameter :: pai = 3.1415926
    real(kind=pr), parameter :: obliquity = 23.43    
    real(kind=pr), parameter :: solar_const = 1365   !w/m2

!------------------------------------------------------------------------------------------------
!start calc total daily average solar radiation (W/m2)
!------------------------------------------------------------------------------------------------
    day_angle  = 2*pai*(day-1)/nday
    dist_ratio = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + &
                    & 0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle) 
    sun_point  = obliquity*cos(2*pai*(173-day)/nday)

    do ilat = 1, nlat, 1
        term = -(tan(pai*sun_point/180.0)*tan(pai*lat(ilat)/180.0))
        where (term<-1)
            term = -1 !polar day
        end where
        where (term>1)
            term = 1  !polar night
        end where
        hour_angle = acos(term)
        all_sw(:,ilat) = (solar_const*dist_ratio/pai)*(hour_angle*sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)+ &
                            & cos(pai*lat(ilat)/180.0)*cos(pai*sun_point/180.0)*sin(hour_angle))
        zenith(:,ilat) = asin(sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)) !used to calc optical depth
    end do

!------------------------------------------------------------------------------------------------
!absorption of shortwave due to water, co2 and o3
!------------------------------------------------------------------------------------------------

end subroutine rad_swf
