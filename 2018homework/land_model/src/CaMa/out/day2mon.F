      program day2mon
! ================================================
      implicit none
! === Basic Parameters ===========================
! Resolution
      integer             ::  nx                  !! number of grids in horizontal
      integer             ::  ny                  !! number of grids in vertical
      integer             ::  nflp
! Input
      character*256       ::  outdir, varname, cyear, buf
      integer             ::  iyear, imon, iday, nday
! Index
      integer             ::  ix, iy
      integer             ::  irec
! Var
      real,allocatable    ::  varbin(:,:), varmon(:,:)
!
      real                ::  rmis
      parameter              (rmis=1.e20)
! File
      character*256       ::  fbin, fmon, fmapdata
      character*256       ::  csufbin, csufmon
      parameter              (csufbin='.bin')
      parameter              (csufmon='.mon')
! Function
      integer             ::  igetday
! ================================================
      call getarg(1,outdir)
      call getarg(2,varname)
      call getarg(3,cyear)
      read(cyear,*) iyear

      fmapdata=trim(outdir)//'mapdata.txt'
      open(11,file=fmapdata,form='formatted')
      read(11,*) buf, nx
      read(11,*) buf, ny
      read(11,*) buf, nflp
      close(11)

      allocate(varbin(nx,ny))
      allocate(varmon(nx,ny))

      fbin=trim(outdir)//trim(varname)//trim(cyear)//trim(csufbin)
      fmon=trim(outdir)//trim(varname)//trim(cyear)//trim(csufmon)
      open(11,file=fbin,form='unformatted',access='direct',recl=4*nx*ny,status='old')
      open(21,file=fmon,form='unformatted',access='direct',recl=4*nx*ny)
      irec=0

      do imon=1, 12
        nday=igetday(iyear,imon)

        varmon(:,:)=0
        do iday=1, nday
          irec=irec+1
          read(11,rec=irec) varbin
          varmon(:,:)=varmon(:,:)+varbin(:,:)
        end do
        varmon(:,:)=varmon(:,:)/real(nday)

        do iy=1, ny
          do ix=1, nx
            if( varbin(ix,iy)==rmis ) varmon(ix,iy)=rmis
          end do
        end do
        write(21,rec=imon) varmon
      end do

      close(11)
      close(12)


      end program day2mon



!!==================================================
      integer function igetday(iyear,imon)
      implicit none
      integer :: iyear, imon

      if( imon==4 .or. imon==6 .or. imon==9 .or. imon==11 )then
        igetday=30
      elseif( imon==2 )then
        igetday=28
        if ( mod(iyear,400)==0 .or. (mod(iyear,100)/=0 .and. mod(iyear,4)==0 )) igetday=29
      else
        igetday=31
      endif

      end function igetday


