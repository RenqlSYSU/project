module create_outcdf_mod
contains 
subroutine create_outcdf(iyyyy)
#ifdef UseCDF
! ===============================================
! create output netcdf file 
!     --------------
!        original : 2011-05-04
!        e. dutra, march 2014, netcdf4
! ===============================================
use parkind1   ,only: jpim, jprb
use mod_input  ,only: lognam,     nx,ny, loutyyyy, coutdir, &
                      coutdir, isyear,ismon,isday,imis,rmis
use mod_output, only: cdfunits,coutvardir,coutvn,coutln,coutunit,&
                      csufcdf, nvars
use mod_map,    only: d1lon,d1lat
use mod_time   ,only: kstep
use netcdf
use lib_netcdf_util

implicit none

integer(kind=jpim)  :: iyyyy

character(len=256)  :: cfile,ctime,cyyyy
integer(kind=jpim)  :: timeid,varid,latid,lonid
integer(kind=jpim)  :: jf


!!===========================

!!=== 
if (kstep > 1) then
  do jf=1, nvars-1
    if ( cdfunits(jf) .eq. 0 ) cycle 
    write(lognam,*) "create_outcdf: closing previous units"
    call ncerror( nf90_close(cdfunits(jf)) )
  enddo
endif

write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a6)') 'seconds since ',isyear,'-',ismon,'-',isday,' 00:00'
write(cyyyy,'(i4)') iyyyy

do jf=1,nvars-1  ! number of output files 

  cfile=trim(coutdir)//trim(coutvardir(jf))//trim(coutvn(jf))//trim(csufcdf)
  if( loutyyyy ) cfile=trim(coutdir)//trim(coutvardir(jf))//trim(coutvn(jf))//trim(cyyyy)//trim(csufcdf)
  if ( trim(coutvardir(jf)) .eq. "NONE" ) cfile="NONE"

  if ( trim(cfile) .eq. "NONE")  then
     cdfunits(jf)=0
  else
    call ncerror( nf90_create(cfile,nf90_netcdf4,cdfunits(jf)),'creating file:'//trim(cfile) )

! === set dimension ===

    call ncerror( nf90_def_dim(cdfunits(jf), 'time', nf90_unlimited, timeid) )

    call ncerror( nf90_def_dim(cdfunits(jf), 'lat', ny, latid) )
    call ncerror( nf90_def_dim(cdfunits(jf), 'lon', nx, lonid) )

! === define variables ===

    call ncerror( nf90_def_var(cdfunits(jf), 'lat', nf90_float, (/latid/), varid) )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'long_name','latitude') )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'units','degrees_north') )

    call ncerror( nf90_def_var(cdfunits(jf), 'lon', nf90_float, (/lonid/), varid) )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'long_name','longitude') )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'units','degrees_east') )

    call ncerror( nf90_def_var(cdfunits(jf), 'time', nf90_double, (/timeid/), varid) ) 
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'long_name','time') )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'units',ctime) )

    
    call ncerror( nf90_def_var(cdfunits(jf), coutvn(jf), nf90_float, (/lonid,latid,timeid/), varid, &
                  deflate_level=6,shuffle=.true.), 'creating variable')
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'long_name', coutln(jf)) )
    call ncerror( nf90_put_att(cdfunits(jf), varid, 'units',     coutunit(jf)) )
    call ncerror( nf90_put_att(cdfunits(jf), varid, '_FillValue',rmis) )

    call ncerror( nf90_enddef(cdfunits(jf)) )

! === put lon lat info ===

    call ncerror ( nf90_inq_varid(cdfunits(jf),'lon',varid),'getting id' )
    call ncerror( nf90_put_var(cdfunits(jf),varid,d1lon))

    call ncerror ( nf90_inq_varid(cdfunits(jf),'lat',varid),'getting id' )
    call ncerror( nf90_put_var(cdfunits(jf),varid,d1lat))


    write(lognam,*) 'cfile: ',trim(cfile),' cvar:',trim(coutvn(jf)),&
                    ' clname: ',trim(coutln(jf)),' cunits: ',trim(coutunit(jf))
    write(lognam,*) 'open in unit: ',cdfunits(jf)
  endif

enddo
#endif
end subroutine create_outcdf

end module create_outcdf_mod
