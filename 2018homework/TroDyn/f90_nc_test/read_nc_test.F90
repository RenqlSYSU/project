program read_nc_test
    implicit none
    include '/usr/netcdf-4331/include/netcdf.inc'

    integer ncid, ierr
    character(len=50) filename

    filename = './CTRL-Clim_month_ave_U.nc'

    PRINT *, filename
    ierr = nf_open(filename, nf_nowrite, ncid)
    print *, ncid
    print *, ierr
    call handle_err(ierr)
end program

    subroutine handle_err(ierr)
        implicit none
        include '/usr/netcdf-4331/include/netcdf.inc'
        integer,intent(in) :: ierr

        if (ierr .NE. nf_noerr) then
            print *, nf_strerror(ierr)
            stop  'stopped'
        endif
    end subroutine handle_err
