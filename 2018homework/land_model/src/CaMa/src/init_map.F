module init_map_mod
contains 
subroutine init_map
! ===============================================
! load map data
! ===============================================
use parkind1   ,only: jpim, jprb, jprm
use mod_input  ,only: nulnam, lognam,tmpnam,nx,ny,nlfp,imis,loutcdf,loutvec,lmapcdf,lmapend,&
                      cnextxy,crivclinc,inpn,pmanriv,pmanfld, &
                      west, east, north, south
use mod_map    ,only: i2nextx,i2nexty,i2rivseq,i2region, regionall, regionthis, &
                      i1seqx,i1seqy,i1next,i2vector,nseqriv,nseqall,nseqmax,rivseqmax, &
                      d1lon, d1lat
use mod_input  ,only: linterp, cinpmat,linterpcdf
use mod_map    ,only: inpx, inpy, inpa
use mod_input  ,only: lpthout,cpthout
use mod_map    ,only: npthout,npthlev,pth_upst,pth_down,pth_dst,pth_elv,pth_wth,pth_man
use mod_output ,only: csufcdf, csufbin, csufvec
#ifdef UseCDF
use netcdf
use lib_netcdf_util
#endif

implicit none

integer(kind=jpim)              :: ix,iy,jx,jy,iseq,ipth,ilev
#ifdef convend
integer(kind=jpim)              :: inpi
#endif
real(kind=jprm),allocatable     :: r2tmp(:,:,:)

#ifdef UseCDF
integer(kind=jpim)              :: ncid,varid
#endif
character*128                   :: cfile1, cfile2

!!=============================
write(lognam,*) '*** init_map: calculate region ***'      

!! *** 1. allocate arrays

allocate( i2nextx(nx,ny) )
allocate( i2nexty(nx,ny) )
allocate( i2rivseq(nx,ny) )
allocate( i2region(nx,ny) )
allocate( i2vector(nx,ny) )
allocate( d1lon(nx) )
allocate( d1lat(ny) )
call read_map

do ix=1,nx
  d1lon(ix)=west+(dble(ix)-0.5)*(east-west)/real(nx)
enddo
do iy=1,ny
  d1lat(iy)=north-(dble(iy)-0.5)*(north-south)/real(ny)
enddo

write(lognam,*) 'init_map: calculate region'      !! decide calculation regions for mpi run
call calc_region

write(lognam,*) 'init_map: calculate 1d river sequence'
allocate(i1seqx(nseqmax))
allocate(i1seqy(nseqmax))
allocate(i1next(nseqmax))

i1seqx(:)=0
i1seqy(:)=0
i1next(:)=0

call calc_1d_seq                                  !! 2d map to 1d vector conversion. for faster calculation

write(lognam,*) 'init_map: nseqriv=',nseqriv
write(lognam,*) 'init_map: nseqall=',nseqall

if( regionthis==1 )then
  write(lognam,*) 'init_map: nseqmax=',nseqmax
  write(lognam,*) 'init_map: rivseqmax=',rivseqmax
endif

if( loutvec )then
  write(lognam,*) 'init_map: wrte_map2vec'
  call wrte_map2vec
endif

!! write map data                                       !! used for combining mpi distributed output into one map
if( regionthis==1 )then
  open(tmpnam,file='./mapdata.txt',form='formatted')
  write(tmpnam,*) 'nx',        nx
  write(tmpnam,*) 'ny',        ny
  write(tmpnam,*) 'nlfp',      nlfp
  write(tmpnam,*) 'regionall', regionall
  write(tmpnam,*) 'nseqmax',   nseqmax
  close(tmpnam)
endif

                     !! runoff interpolation
if( linterp )then
  write(lognam,*) 'init_map: input matrix ', cinpmat
  write(lognam,*) '  nx, ny, inpn =', nx, ny, inpn
  allocate(inpx(nx,ny,inpn),inpy(nx,ny,inpn),inpa(nx,ny,inpn),r2tmp(nx,ny,inpn))
  if ( linterpcdf ) then 
#ifdef UseCDF
    call ncerror (nf90_open(cinpmat,nf90_nowrite,ncid),'opening '//trim(cinpmat) )

    write(lognam,*)'init_map: area:',trim(cinpmat)
    call ncerror ( nf90_inq_varid(ncid,'area',varid),'getting id' )
    call ncerror ( nf90_get_var(ncid,varid,inpa,(/1,1,1/),(/nx,ny,inpn/)),'reading data' ) 

    write(lognam,*)'init_map: indx:',trim(cinpmat)
    call ncerror ( nf90_inq_varid(ncid,'indx',varid),'getting id' )
    call ncerror ( nf90_get_var(ncid,varid,inpx,(/1,1,1/),(/nx,ny,inpn/)),'reading data' ) 

    write(lognam,*)'init_map: indy:',trim(cinpmat)
    call ncerror ( nf90_inq_varid(ncid,'indy',varid),'getting id' )
    call ncerror ( nf90_get_var(ncid,varid,inpy,(/1,1,1/),(/nx,ny,inpn/)),'reading data' ) 
#endif
  else
    open(tmpnam,file=cinpmat,form='unformatted',access='direct',recl=4*nx*ny*inpn)
    read(tmpnam,rec=1) inpx
    read(tmpnam,rec=2) inpy
    read(tmpnam,rec=3) r2tmp

#ifdef convend
    if ( lmapend )then
      do iy=1, ny
        do ix=1, nx
          do inpi=1, inpn
            call endian4(inpx(ix,iy,inpi))
            call endian4(inpx(ix,iy,inpi))
            call endian4(r2tmp(ix,iy,inpi))
          end do
        end do
      end do
    endif
#endif

    inpa(:,:,:)=dble(r2tmp(:,:,:))
    close(tmpnam)
    deallocate(r2tmp)
  endif
endif


                     !! bifurcation channel
if( lpthout )then
  open(tmpnam,file=cpthout,form='formatted')
  read(tmpnam,*) npthout,npthlev

  allocate(pth_upst(npthout))
  allocate(pth_down(npthout))

  allocate(pth_dst(npthout))
  allocate(pth_elv(npthout,npthlev))
  allocate(pth_wth(npthout,npthlev))
  allocate(pth_man(npthlev))

  do ipth=1, npthout
    read(tmpnam,*) ix, iy, jx, jy, pth_dst(ipth), ( pth_elv(ipth,ilev),pth_wth(ipth,ilev),ilev=1,npthlev )
    pth_upst(ipth)=i2vector(ix,iy)
    pth_down(ipth)=i2vector(jx,jy)
  end do

  do ilev=1, npthlev
    if( ilev==1 )then
      pth_man(ilev)=pmanriv
    else
      pth_man(ilev)=pmanfld
    endif
  end do

  close(tmpnam)
endif

contains


! ================================================
subroutine read_map

if ( .not. lmapcdf ) then
  write(lognam,*)'init_map: nextxy file: ',trim(cnextxy)
  open(tmpnam,file=cnextxy,form='unformatted',access='direct',recl=4*nx*ny)
  read(tmpnam,rec=1) i2nextx
  read(tmpnam,rec=2) i2nexty
  close(tmpnam)

#ifdef convend
  if ( lmapend )then
    do iy=1, ny
      do ix=1, nx
        call endian4(i2nextx(ix,iy))
        call endian4(i2nexty(ix,iy))
      end do
    end do
  endif
#endif

else
#ifdef UseCDF
  call ncerror (nf90_open(crivclinc,nf90_nowrite,ncid),'opening '//trim(crivclinc) )

  write(lognam,*)'init_map: nextx:',trim(crivclinc)
  call ncerror ( nf90_inq_varid(ncid,'nextx',varid),'getting id' )
  call ncerror ( nf90_get_var(ncid,varid,i2nextx),'reading data' ) 

  write(lognam,*)'init_map: nexty:',trim(crivclinc)
  call ncerror ( nf90_inq_varid(ncid,'nexty',varid),'getting id' )
  call ncerror ( nf90_get_var(ncid,varid,i2nexty),'reading data' )

  write(lognam,*)'init_map: lat:',trim(crivclinc)
  call ncerror ( nf90_inq_varid(ncid,'lat',varid),'getting id' )
  call ncerror ( nf90_get_var(ncid,varid,d1lat),'reading data' )

  write(lognam,*)'init_map: lon:',trim(crivclinc)
  call ncerror ( nf90_inq_varid(ncid,'lon',varid),'getting id' )
  call ncerror ( nf90_get_var(ncid,varid,d1lon),'reading data' )

  call ncerror( nf90_close(ncid))
#endif
endif

call calc_rivseq

end subroutine read_map


! ================================================
subroutine calc_rivseq
implicit none
integer(kind=jpim) :: ix, iy, jx, jy, seqnow, next

write(lognam,*) '  calculate river seqence'  
i2rivseq(:,:)=1
do iy=1, ny
  do ix=1, nx
    if( i2nextx(ix,iy)>0 )then
      jx=i2nextx(ix,iy)
      jy=i2nexty(ix,iy)
      i2rivseq(jx,jy)=0
    elseif( i2nextx(ix,iy)==-9999 )then
      i2rivseq(ix,iy)=-9999
    endif
  end do
end do

seqnow=0
next=1
do while(next>0 )
  seqnow=seqnow+1
  next=0
  do iy=1, ny
    do ix=1, nx
      if( i2rivseq(ix,iy)==seqnow )then
        if( i2nextx(ix,iy)>0 )then
          jx=i2nextx(ix,iy)
          jy=i2nexty(ix,iy)
          if( i2rivseq(jx,jy)<=seqnow )then
            next=next+1
            i2rivseq(jx,jy)=seqnow+1
          endif
        endif
      endif
    end do
  enddo
end do
rivseqmax=seqnow
write(lognam,*) '  rivseqmax = ', rivseqmax


end subroutine calc_rivseq





!!==================================================
subroutine calc_region    !! evenly allocate pixels to mpi nodes
implicit none

integer(kind=jpim),allocatable :: i2basin(:,:)
integer(kind=jpim),allocatable :: i2upgrid(:,:)

integer(kind=jpim),allocatable :: basingrid(:)
integer(kind=jpim),allocatable :: regiongrid(:)
integer(kind=jpim),allocatable :: basinregion(:)

integer(kind=jpim) :: ix,iy,jx,jy,kx,ky,upg, basinmax, iregion, jregion, gridmin, ibasin

!!==============================

allocate(i2basin(nx,ny))
allocate(i2upgrid(nx,ny))

write(lognam,*)'  calculate upstream grid number'
i2upgrid(:,:)=0
do iy=1, ny
  do ix=1, nx
    if( i2rivseq(ix,iy)==1 ) then
      jx=ix
      jy=iy
      i2upgrid(jx,jy)=1
      upg=i2upgrid(jx,jy)
      do while( i2nextx(jx,jy)>0 )  !! if river reaches mouth, end loop
        kx=i2nextx(jx,jy)
        ky=i2nexty(jx,jy)
        jx=kx
        jy=ky
        if( i2upgrid(jx,jy)==0 )then               !! grids firstly checked
          i2upgrid(jx,jy)=upg+1
          upg=i2upgrid(jx,jy)
        else                                       !! grids already checked
          i2upgrid(jx,jy)=i2upgrid(jx,jy)+upg
        endif
      end do
    endif
  end do
end do

write(lognam,*)'  calculate basin'
i2basin(:,:)=0
ibasin=0
do iseq=rivseqmax, 1, -1
  do iy=1, ny
    do ix=1, nx
      if( i2rivseq(ix,iy)==iseq .and. i2nextx(ix,iy)<0 .and. i2rivseq(ix,iy)/=imis )then
        ibasin=ibasin+1
        i2basin(ix,iy)=ibasin
      endif
    end do
  end do
end do
basinmax=ibasin

do iy=1, ny
  do ix=1, nx
    jx=ix
    jy=iy
    do while( i2basin(jx,jy)==0 .and. i2nextx(jx,jy)>0 )
      kx=i2nextx(jx,jy)
      ky=i2nexty(jx,jy)
      jx=kx
      jy=ky
    end do
    ibasin=i2basin(jx,jy)
    jx=ix
    jy=iy
    do while( i2basin(jx,jy)==0 .and. i2nextx(jx,jy)>0  )
      i2basin(jx,jy)=ibasin
      kx=i2nextx(jx,jy)
      ky=i2nexty(jx,jy)
      jx=kx
      jy=ky
    end do
  end do
end do

allocate(basingrid(basinmax))
allocate(basinregion(basinmax))
allocate(regiongrid(regionall))

write(lognam,*)'  allocate basin to cpu (mpi)'
do iy=1, ny
  do ix=1, nx
    if( i2nextx(ix,iy)<0 .and. i2nextx(ix,iy)/=imis )then
      ibasin=i2basin(ix,iy)
      basingrid(ibasin)=i2upgrid(ix,iy)
    endif
  end do
end do


regiongrid(:)=0
jregion=1
do ibasin=1, basinmax
  gridmin=nx*ny
  do iregion=1, regionall
    if( regiongrid(iregion) < gridmin )then
      gridmin=regiongrid(iregion)
      jregion=iregion
    endif
  end do
  basinregion(ibasin)=jregion
  regiongrid(jregion)=regiongrid(jregion)+basingrid(ibasin)
end do

nseqmax=0
do iregion=1, regionall
  if( regionthis==1 )then
    write(nulnam,*), 'calc_region: ', iregion, regiongrid(iregion)
  endif
  nseqmax=max(nseqmax,regiongrid(iregion))
end do

do iy=1, ny
  do ix=1, nx
    if( i2basin(ix,iy)>0 .and. i2nextx(ix,iy)/=-9999 )then
      ibasin=i2basin(ix,iy)
      i2region(ix,iy)=basinregion(ibasin)
    endif
  end do
end do

deallocate(basingrid)
deallocate(basinregion)
deallocate(regiongrid)

deallocate(i2basin)
deallocate(i2upgrid)

end subroutine calc_region


! ================================================
subroutine calc_1d_seq
implicit none

i2vector(:,:)=0
iseq=0
do iy=1, ny
  do ix=1, nx
    if( i2nextx(ix,iy).gt.0 .and. i2region(ix,iy)==regionthis )then
      iseq=iseq+1
      i1seqx(iseq)=ix
      i1seqy(iseq)=iy
      i2vector(ix,iy)=iseq
    endif
  end do
end do
nseqriv=iseq

do iy=1, ny
  do ix=1, nx
    if( i2nextx(ix,iy).lt.0 .and. i2nextx(ix,iy).ne.imis .and. i2region(ix,iy)==regionthis )then
      iseq=iseq+1
      i1seqx(iseq)=ix
      i1seqy(iseq)=iy
      i2vector(ix,iy)=iseq
    endif
  end do
end do
nseqall=iseq

do iseq=1, nseqall
  ix=i1seqx(iseq)
  iy=i1seqy(iseq)
  if( i2nextx(ix,iy)>0 )then
    jx=i2nextx(ix,iy)
    jy=i2nexty(ix,iy)
    i1next(iseq)=i2vector(jx,jy)
  else
    i1next(iseq)=-9
  endif
end do

end subroutine calc_1d_seq


!!==================================================
subroutine wrte_map2vec       !! 1d sequence vector informtion required to convert vector output to 2d map

implicit none
#ifdef UseCDF
integer(kind=jpim) :: xid,yid,landid
#endif

if( loutcdf )then

  cfile1='./region'//trim(csufcdf)
  cfile2='./ind_xy'//trim(csufcdf)

#ifdef UseCDF
  if( regionthis==1 )then

    call ncerror( nf90_create(cfile1, nf90_64bit_offset, ncid) )
    !!== dimensions
    call ncerror( nf90_def_dim(ncid,'lon',   nx,   xid) )
    call ncerror( nf90_def_dim(ncid,'lat',   ny,   yid) )
    !!== variables
    call ncerror( nf90_def_var(ncid, 'lon', nf90_float, (/xid/), varid) ) !done
    call ncerror( nf90_put_att(ncid, varid, 'long_name','longitude') )
    call ncerror( nf90_put_att(ncid, varid, 'units','degrees_east') )

    call ncerror( nf90_def_var(ncid, 'lat', nf90_float, (/yid/), varid) ) !done
    call ncerror( nf90_put_att(ncid, varid, 'long_name','latitude') )
    call ncerror( nf90_put_att(ncid, varid, 'units','degrees_north') )

    call ncerror( nf90_def_var(ncid, 'region', nf90_int, (/xid,yid/), varid) ) ! done
    call ncerror( nf90_put_att(ncid, varid, 'long_name','mpi region') )
    call ncerror( nf90_put_att(ncid, varid, 'units','region index') )

    call ncerror( nf90_enddef(ncid) )
    !!== put vars 
    call ncerror( nf90_inq_varid(ncid,'lon',varid))
    call ncerror( nf90_put_var(ncid,varid,d1lon))

    call ncerror( nf90_inq_varid(ncid,'lat',varid))
    call ncerror( nf90_put_var(ncid,varid,d1lat))

    call ncerror( nf90_inq_varid(ncid,'region',varid))
    call ncerror( nf90_put_var(ncid,varid,i2region))

    call ncerror( nf90_close(ncid) )
  endif

  call ncerror( nf90_create(cfile2,nf90_64bit_offset,ncid) )
  !!== dimensions
  call ncerror( nf90_def_dim(ncid, 'land', nseqmax, landid) )
  call ncerror( nf90_def_dim(ncid,'lon',   nx,   xid) )
  call ncerror( nf90_def_dim(ncid,'lat',   ny,   yid) )
    !!== variables
  call ncerror( nf90_def_var(ncid, 'lon', nf90_float, (/xid/), varid) ) !done
  call ncerror( nf90_put_att(ncid, varid, 'long_name','longitude') )
  call ncerror( nf90_put_att(ncid, varid, 'units','degrees_east') )

  call ncerror( nf90_def_var(ncid, 'lat', nf90_float, (/yid/), varid) ) !done
  call ncerror( nf90_put_att(ncid, varid, 'long_name','latitude') )
  call ncerror( nf90_put_att(ncid, varid, 'units','degrees_north') )

  call ncerror( nf90_def_var(ncid, 'land', nf90_int, (/landid/), varid) ) ! done
  call ncerror( nf90_put_att(ncid, varid, 'compress','y x') )
  call ncerror( nf90_put_att(ncid, varid, 'info','y/x index 1:ny(nx=1),ny+1:ny*2(nx=2),etc') )
  
  call ncerror( nf90_def_var(ncid, 'land_ind', nf90_int, (/xid,yid/), varid) ) ! done
  call ncerror( nf90_put_att(ncid, varid, 'long_name','land index for x,y pairs') )
   
  call ncerror( nf90_def_var(ncid, 'ind_x', nf90_int, (/landid/), varid) ) !done
  call ncerror( nf90_put_att(ncid, varid, 'compress','index x ') )

  call ncerror( nf90_def_var(ncid, 'ind_y', nf90_int, (/landid/), varid) ) !done
  call ncerror( nf90_put_att(ncid, varid, 'compress','index y ') )

  call ncerror( nf90_enddef(ncid) )
  !!== put vars
  call ncerror( nf90_inq_varid(ncid,'ind_x',varid))
  call ncerror( nf90_put_var(ncid,varid,i1seqx))

  call ncerror( nf90_inq_varid(ncid,'ind_y',varid))
  call ncerror( nf90_put_var(ncid,varid,i1seqy))
  
  call ncerror( nf90_inq_varid(ncid,'land_ind',varid))
  call ncerror( nf90_put_var(ncid,varid,i2vector))

  call ncerror( nf90_close(ncid) )
#endif

else !! binary output
  cfile1='./region'//trim(csufbin)
  cfile2='./ind_xy'//trim(csufvec)

  if( regionthis==1 )then
    open(tmpnam,file=cfile1,form='unformatted',access='direct',recl=4*nx*ny)
    write(tmpnam,rec=1) i2region
    close(tmpnam)
  endif
  open(tmpnam,file=cfile2,form='unformatted',access='direct',recl=4*nseqmax)
  write(tmpnam,rec=1) i1seqx
  write(tmpnam,rec=2) i1seqy
  close(tmpnam)
endif

end subroutine wrte_map2vec



end subroutine init_map
end module init_map_mod
