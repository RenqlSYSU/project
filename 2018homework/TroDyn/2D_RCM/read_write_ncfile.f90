!this module contains the subroutines used to write and read nc file 
!by renql 20190412

module write_read_ncfile
    implicit none
    include '/usr/netcdf-4331/include/netcdf.inc'

contains

    subroutine write_nc(filename,day,lev,lat,var)
        implicit none
    end subroutine write_nc

    subroutine read_nc(filename,day,lev,lat,var)
        implicit none
    end subroutine read_nc

!-----used to show error status-----------------------------
    subroutine handle_err(status)
        implicit none
        integer  status

        if (status .NE. nf_noerr) then
            print *, nf_strerror(status)
            stop  'stopped'
        endif
    end subroutine handle_err

end module write_read_ncfile
