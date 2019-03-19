      program txt_vector
! ================================================
      implicit none
! Resolution
      integer             ::  ix, iy, jx, jy
      integer             ::  nx, ny                        !! grid numbers
      real                ::  west0, north0                 !! map west and north edge
      real                ::  gsize                         !! map grid size
! domain
      real                ::  west, east, north, south      !! plot domain
! vars
      integer,allocatable ::  nextx(:,:), nexty(:,:)        !! downstream xy
      real,allocatable    ::  lon(:,:), lat(:,:)
      real,allocatable    ::  uparea(:,:)                   !! dorainage area [m2]
      real                ::  lon1, lon2, lat1, lat2
! files
      character*128       ::  params
      character*128       ::  rfile1, rfile2, rfile3
      parameter              (params='./map/params.txt')
      parameter              (rfile1='./map/nextxy.bin')
      parameter              (rfile2='./map/lonlat.bin')
      parameter              (rfile3='./map/uparea.bin')

      character*64        ::  buf
! ================================================
      call getarg(1,buf)
      read(buf,*) west
      call getarg(2,buf)
      read(buf,*) east
      call getarg(3,buf)
      read(buf,*) north
      call getarg(4,buf)
      read(buf,*) south

      open(10,file=params,form='formatted')
      read(10,*) west0
      read(10,*) north0
      read(10,*) nx
      read(10,*) ny
      read(10,*) gsize
      close(10)

      allocate(nextx(nx,ny),nexty(nx,ny))
      allocate(lon(nx,ny),lat(nx,ny))
      allocate(uparea(nx,ny))

      open(11,file=rfile1,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) nextx
      read(11,rec=2) nexty
      close(11)

      open(11,file=rfile2,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) lon
      read(11,rec=2) lat
      close(11)

      open(11,file=rfile3,form='unformatted',access='direct',recl=4*nx*ny)
      read(11,rec=1) uparea
      close(11)

      do iy=1, ny
        do ix=1, nx
          if( nextx(ix,iy)>0 )then
            jx=nextx(ix,iy)
            jy=nexty(ix,iy)
            lon1=lon(ix,iy)
            lon2=lon(jx,jy)
            lat1=lat(ix,iy)
            lat2=lat(jx,jy)

            if( lon1>=west .and. lon1<=east .and. lat1<=north .and. lat1>=south )then
              print '(4f12.5,f12.1)', lon1, lat1, lon2, lat2, uparea(ix,iy)/1000.**2
            elseif( lon2>=west .and. lon2<=east .and. lat2<=north .and. lat2>=south )then
              print '(4f12.5,f12.1)', lon2, lat2, lon1, lat1, uparea(ix,iy)/1000.**2
            endif
          elseif( nextx(ix,iy)==-9 .or. nextx(ix,iy)==-10 )then
            print '(4f12.5,f12.1)', lon(ix,iy), lat(ix,iy), -999.,      -999.,      &
                                                   uparea(ix,iy)/1000.**2
          endif
        end do
      end do
      
! === Finish =====================================
      end program txt_vector
