module control_rest_mod
contains 
subroutine control_rest(iyyyy,imm,idd)
! ===============================================
! control restart files
! ===============================================
use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: lognam, tmpnam, nx,ny, crestdir, lrestcdf, lstoonly, restfreq,&
                      isyear, ismon,  isday, rmis, dmis,lpthout
use mod_time   ,only: kstep, nsteps, kmin, kminstart, khhmm
use mod_output ,only: cvnrestart,csufbin,csufcdf
use mod_map    ,only: nseqmax, npthout, npthlev,d1lon,d1lat
use mod_prog   ,only: d2rivsto_out, d2fldsto_out, d2rivout_pre, d2fldout_pre, d1pthflw_pre, d2rivdph_pre, d2fldsto_pre

implicit none

!!== input
integer(kind=jpim),intent(in)  :: iyyyy,imm,idd

! ===============================================
if ( .not. lrestcdf )then
  if ( restfreq==0 .and. imm==1 .and. idd == 1 .and. khhmm == 0 ) call wrte_rest_bin  ! end of year
  if ( restfreq==1 .and. khhmm .eq. 0 )                           call wrte_rest_bin  ! daily
  if ( restfreq==2 .and. idd == 1 .and. khhmm == 0 )              call wrte_rest_bin  ! end of month
elseif ( lrestcdf ) then ! netcdf restart 
  !if ( restfreq .eq. 0 .and. kstep .eq. nsteps ) call wrte_rest_cdf  ! final of run
  if ( restfreq .eq. 0 .and. imm .eq. 1 .and. idd .eq. 1 .and. khhmm .eq. 0 ) call wrte_rest_bin  ! end of year
  if ( restfreq .eq. 1 .and. khhmm .eq. 0 )                                   call wrte_rest_cdf  ! daily 
  if ( restfreq .eq. 2 .and. idd .eq. 1 .and. khhmm .eq. 0 )                  call wrte_rest_cdf  ! end of month 
end if 

contains





! ==================================================
subroutine wrte_rest_bin
implicit none
! 
real(kind=jprm)       :: r2temp(nx,ny)
real(kind=jprm)       :: r1pth(npthout,npthlev)
character(len=256)    :: cfile,cdate


!!==============================

write(cdate,'(i4.4,2i2.2)') iyyyy, imm, idd

cfile=trim(crestdir)//trim(cvnrestart)//trim(cdate)//trim(csufbin)
write(lognam,*) 'wrte_rest: write restart bin:',cfile

open(tmpnam,file=cfile,form='unformatted',access='direct',recl=4*nx*ny)
  call vec2map(d2rivsto_out,r2temp)
   write(tmpnam,rec=1) r2temp
  call vec2map(d2fldsto_out,r2temp)
   write(tmpnam,rec=2) r2temp
  if( .not. lstoonly )then
    call vec2map(d2rivout_pre,r2temp)
     write(tmpnam,rec=3) r2temp
    call vec2map(d2fldout_pre,r2temp)
     write(tmpnam,rec=4) r2temp
    call vec2map(d2rivdph_pre,r2temp)
     write(tmpnam,rec=5) r2temp
    call vec2map(d2fldsto_pre,r2temp)
     write(tmpnam,rec=6) r2temp
  endif
close(tmpnam)

if( lpthout )then
  if( .not. lstoonly )then
    cfile=trim(crestdir)//trim(cvnrestart)//trim(cdate)//trim(csufbin)//'.pth'
    write(lognam,*) 'wrte_rest: write restart bin:',cfile

    open(tmpnam,file=cfile,form='unformatted',access='direct',recl=4*npthout*npthlev)
    r1pth(:,:)=real(d1pthflw_pre(:,:))
    write(tmpnam,rec=1) r1pth
    close(tmpnam)
  endif
endif

end subroutine wrte_rest_bin







! ==================================================
subroutine wrte_rest_cdf
#ifdef UseCDF
use netcdf
use lib_netcdf_util
implicit none
! 
character(len=256)             :: cfile,cdate,ctime,cvar
integer(kind=jpim)             :: ncid,varid,latid,lonid,timeid,jf
real(kind=jprb)                :: xtime ! seconds since start of the run ! 
real(kind=jprb)                :: r2temp(nx,ny)

if( lpthout )then
  write(lognam,*) 'wrte_rest: lpthout is active, but we don not save it to the restart file!'
endif 

xtime=real( (kmin-kminstart),jprb) *60._jprb
write(ctime,'(a14,i4.4,a1,i2.2,a1,i2.2,a6)') 'seconds since ',isyear,'-',ismon,'-',isday,' 00:00'
write(cdate,'(i4.4,i2.2,i2.2)') iyyyy,imm,idd
cfile=trim(crestdir)//trim(cvnrestart)//trim(cdate)//trim(csufcdf)
write(lognam,*) 'wrte_rest:create restart netcdf:',cfile
!! create netcdf file
call ncerror( nf90_create(cfile,nf90_netcdf4,ncid),'creating file:'//trim(cfile) )

!! dimensions 
call ncerror( nf90_def_dim(ncid, 'time', nf90_unlimited, timeid) )
call ncerror( nf90_def_dim(ncid, 'lat', ny, latid) )
call ncerror( nf90_def_dim(ncid, 'lon', nx, lonid) )

!! variables: 
call ncerror( nf90_def_var(ncid, 'lat', nf90_float, (/latid/), varid) )
call ncerror( nf90_put_att(ncid, varid, 'long_name','latitude') )
call ncerror( nf90_put_att(ncid, varid, 'units','degrees_north') )

call ncerror( nf90_def_var(ncid, 'lon', nf90_float, (/lonid/), varid) )
call ncerror( nf90_put_att(ncid, varid, 'long_name','longitude') )
call ncerror( nf90_put_att(ncid, varid, 'units','degrees_east') )

call ncerror( nf90_def_var(ncid, 'time', nf90_double, (/timeid/), varid) ) 
call ncerror( nf90_put_att(ncid, varid, 'long_name','time') )
call ncerror( nf90_put_att(ncid, varid, 'units',ctime) )

call ncerror( nf90_def_var(ncid, 'rivsto', nf90_double, (/lonid,latid,timeid/), varid, &
                          deflate_level=6,shuffle=.true.), 'creating variable')
call ncerror( nf90_put_att(ncid, varid, 'long_name',"river storage" ) )
call ncerror( nf90_put_att(ncid, varid, 'units',"m3") )
call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis),'in here?' )

call ncerror( nf90_def_var(ncid, 'fldsto', nf90_double, (/lonid,latid,timeid/), varid, &
                          deflate_level=6,shuffle=.true.), 'creating variable')  
call ncerror( nf90_put_att(ncid, varid, 'long_name',"flood plain storage" ) )
call ncerror( nf90_put_att(ncid, varid, 'units',"m3") )
call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis) )

if( .not. lstoonly )then

  call ncerror( nf90_def_var(ncid, 'rivout', nf90_double, (/lonid,latid,timeid/), varid, &
                            deflate_level=6,shuffle=.true.), 'creating variable')  
  call ncerror( nf90_put_att(ncid, varid, 'long_name',"river outflow" ) )
  call ncerror( nf90_put_att(ncid, varid, 'units',"m3/s") )
  call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis) )

  call ncerror( nf90_def_var(ncid, 'fldout', nf90_double, (/lonid,latid,timeid/), varid, &
                            deflate_level=6,shuffle=.true.), 'creating variable')  
  call ncerror( nf90_put_att(ncid, varid, 'long_name',"floodplain outflow" ) )
  call ncerror( nf90_put_att(ncid, varid, 'units',"m3/s") )
  call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis) )

  call ncerror( nf90_def_var(ncid, 'rivdph', nf90_double, (/lonid,latid,timeid/), varid, &
                            deflate_level=6,shuffle=.true.), 'creating variable')  
  call ncerror( nf90_put_att(ncid, varid, 'long_name',"river depth" ) )
  call ncerror( nf90_put_att(ncid, varid, 'units',"m") )
  call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis) )

  call ncerror( nf90_def_var(ncid, 'fldsto_pre', nf90_double, (/lonid,latid,timeid/), varid, &
                            deflate_level=6,shuffle=.true.), 'creating variable')  
  call ncerror( nf90_put_att(ncid, varid, 'long_name',"floodplain storage prev" ) )
  call ncerror( nf90_put_att(ncid, varid, 'units',"m3") )
  call ncerror( nf90_put_att(ncid, varid, '_FillValue',dmis) )

! call ncerror( nf90_def_var(ncid, 'pthflw', nf90_double, (/landid,timeid/), varid) ) 
! call ncerror( nf90_put_att(ncid, varid, 'long_name',"floodpath outflow" ) )
! call ncerror( nf90_put_att(ncid, varid, 'units',"m3/s") )

endif

call ncerror( nf90_enddef(ncid) )


call ncerror( nf90_inq_varid(ncid,'time',varid))
call ncerror( nf90_put_var(ncid,varid,xtime) )

call ncerror ( nf90_inq_varid(ncid,'lon',varid),'getting id' )
call ncerror( nf90_put_var(ncid,varid,d1lon))

call ncerror ( nf90_inq_varid(ncid,'lat',varid),'getting id' )
call ncerror( nf90_put_var(ncid,varid,d1lat))

do jf=1,6
  select case(jf)
    case (1)
        cvar='rivsto'
        call vec2mapd(d2rivsto_out,r2temp)
    case (2)
        cvar='fldsto'
        call vec2mapd(d2fldsto_out,r2temp)
    case (3)
        cvar='rivout'
        call vec2mapd(d2rivout_pre,r2temp)
    case (4)
        cvar='fldout'
        call vec2mapd(d2fldout_pre,r2temp)
    case (5)
        cvar='rivdph'
        call vec2mapd(d2rivdph_pre,r2temp)
    case (6)
        cvar='fldsto_pre'
        call vec2mapd(d2fldsto_pre,r2temp)
  end select

  if( (.not. lstoonly) .or. (jf<=2) )then
    call ncerror( nf90_inq_varid(ncid,trim(cvar),varid))
    call ncerror( nf90_put_var(ncid,varid,r2temp,(/1,1,1/),(/nx,ny,1/)) )
  endif
enddo

call ncerror( nf90_sync(ncid) )
call ncerror( nf90_close(ncid ) )

 
write(lognam,*) 'wrte_rest: write restart netcdf:',cfile
 
#endif
end subroutine wrte_rest_cdf






! ================================================
subroutine vec2map(d2vec,r2map)
use mod_map    ,only: i2vector
implicit none

real(kind=jprb)              :: d2vec(nseqmax,1)
real(kind=jprm)              :: r2map(nx,ny)

!$ save
integer(kind=jpim)           ::  ix,iy,iseq
!$omp threadprivate(ix,iseq)
! ======

!$omp parallel do
do iy=1, ny
  do ix=1, nx
    if (i2vector(ix,iy)>0) then
      iseq=i2vector(ix,iy)
      r2map(ix,iy) = real(d2vec(iseq,1))
    else
      r2map(ix,iy) = rmis
    end if
  end do
end do
!$omp end parallel do

end subroutine vec2map

! ================================================
subroutine vec2mapd(d2vec,r2map)
use mod_map    ,only: i2vector
implicit none

real(kind=jprb)              :: d2vec(nseqmax,1)
real(kind=jprb)              :: r2map(nx,ny)

!$ save
integer(kind=jpim)           ::  ix,iy,iseq
!$omp threadprivate(ix,iseq)
! ======
r2map(:,:) = dmis
!$omp parallel do
do iy=1, ny
  do ix=1, nx
    if (i2vector(ix,iy)>0) then
      iseq=i2vector(ix,iy)
      r2map(ix,iy) = d2vec(iseq,1)
    end if
  end do
end do
!$omp end parallel do

end subroutine vec2mapd

end subroutine control_rest
end module control_rest_mod
