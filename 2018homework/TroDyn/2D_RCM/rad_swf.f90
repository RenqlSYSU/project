!This module 
!by renql 20190412
!

subroutine rad_swf(net_flux)
    use constant
    implicit none
    
    integer :: i, ilat, iday !used to assignment or cycling
    real    :: f

    real(kind=pr), dimension(nbndsw) :: sw_irrad !total daily average solar radiation (W/m2)
    real(kind=pr), dimension(nday,nlat,nbndsw) :: toa_sw_irrad !total daily average solar radiation (W/m2)
    real(kind=pr), dimension(nday,nlat) :: zenith !mean zenith angle 
    real(kind=pr), dimension(nday,nlat) :: ratio  !ratio between solar_constant and earth recived

    real(kind=pr), parameter :: obliquity = 23.43    
    real(kind=pr), parameter :: solar_const = 1368  !w/m2
    real(kind=pr), parameter :: solar_temp  = 400  !K, Calculated from the solar constant
    
    real(kind=pr):: day_angle, dist_ratio, sun_point, hour_angle, term

!------------------------------------------------------------------------------------------------
!start calc total daily average solar radiation (W/m2)
!------------------------------------------------------------------------------------------------
    call calc_bnd_irradiance(solar_temp,nbndsw,sw_lght,sw_irrad)

    do iday = 1, nday, 1
        day_angle  = 2*pai*(iday-1)/nday
        dist_ratio = 1.00011 + 0.034221*cos(day_angle) + 0.00128*sin(day_angle) + &
                        & 0.000719*cos(2*day_angle) + 0.000077*sin(2*day_angle) 
        sun_point  = obliquity*cos(2*pai*(173-iday)/nday)
    
        do ilat = 1, nlat, 1
            term = -(tan(pai*sun_point/180.0)*tan(pai*lat(ilat)/180.0))
            where (term<-1)
                term = -1 !polar day
            end where
            where (term>1)
                term = 1  !polar night
            end where
            hour_angle = acos(term)
            toa_sw_irrad(iday,ilat,:) = (sw_irrad*dist_ratio/pai)*(hour_angle*sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)+ &
                                            & cos(pai*lat(ilat)/180.0)*cos(pai*sun_point/180.0)*sin(hour_angle))
            zenith(iday,ilat) = asin(sin(pai*lat(ilat)/180.0)*sin(pai*sun_point/180.0)) !used to calc optical depth
        end do
    end do

!------------------------------------------------------------------------------------------------
!absorption of shortwave due to water, co2 and o3
!------------------------------------------------------------------------------------------------

end subroutine rad_swf
