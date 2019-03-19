      program proj_flood
! ===============================================
      implicit none
! CaMa-Flood parameters       
      character*64             ::  param                         !! river map parameters
      integer                  ::  iXX, iYY
      integer                  ::  nXX, nYY                      !! grid number (river network map)
      integer                  ::  nflp                          !! floodplain layers
      real                     ::  gsize                         !! grid size [deg]
      real                     ::  west, east, north, south      !! domain (river network map)

! HydroSHEDS parameters                                        !! (from location.txt)
      integer                  ::  i, iarea, narea               !! area ID
      character*3              ::  area                          !! area code
      integer                  ::  ix, iy                        
      integer                  ::  nx, ny                        !! grid number (hires data)
      real                     ::  csize                         !! size of pixel [deg]
      real                     ::  lon_ori                       !! west  edge
      real                     ::  lat_ori                       !! north edge

      character*64             ::  list_loc
      character*3,allocatable  ::  list_area(:)                  !! area code
      integer,allocatable      ::  list_nx(:),  list_ny(:)
      real,allocatable         ::  list_lon(:), list_lat(:)
!
      integer*2,allocatable    ::  nextx(:,:)                    !! downstream (jXX,jYY)
      integer*2,allocatable    ::  catmx(:,:), catmy(:,:)        !! catchment (iXX,iYY) of pixel (ix,iy)
      real,allocatable         ::  flddif(:,:)                   !! height above channel [m]

      real,allocatable         ::  flddph(:,:)                   !! flood depth [m] (coarse resolution)
!
      real,allocatable         ::  flood(:,:)                    !! downscaled flood depth [m]
!
      character*64             ::  mapdir, hires
      parameter                   (mapdir='./map/')              !! map directory (please make a symbolic link)
      parameter                   (hires='./map/hires/')         !! output data directory (please make a symbolic link)
      character*128            ::  fnextxy, fgridxy, fflddif, fflood
      character*128            ::  fflddph
      integer                  ::  trec
      character*128            ::  buf
! ===============================================
      call getarg(1,area)       !! downscale area
      call getarg(2,fflddph)    !! downscale file
      call getarg(3,buf)        !! downscale file record number
      if( buf/='' )then
        read(buf,*) trec
      else
        trec=1
      endif


      param=trim(mapdir)//'params.txt'

      open(11,file=param,form='formatted')
      read(11,*) west
      read(11,*) north
      read(11,*) nXX
      read(11,*) nYY
      read(11,*) gsize
      read(11,*) nflp
      read(11,*) narea
      read(11,*) csize
      close(11)

      east =west +real(nXX)*gsize
      south=north-real(nXX)*gsize

! ==========

      list_loc=trim(hires)//'location.txt'
      allocate(list_area(narea),list_nx(narea),list_ny(narea),list_lon(narea),list_lat(narea))

      open(11,file=list_loc,form='formatted')
      read(11,*)
      read(11,*) buf, (list_area(i) ,i=1,narea)
      read(11,*) buf, (list_lon(i)  ,i=1,narea)
      read(11,*) buf, (list_lat(i)  ,i=1,narea)
      read(11,*) buf, (list_nx(i)   ,i=1,narea)
      read(11,*) buf, (list_ny(i)   ,i=1,narea)
      close(11)

      nx=-9999
      do i=1, narea
        if( area==list_area(i) )then
          iarea=i
          nx=list_nx(i)
          ny=list_ny(i)
          lon_ori=list_lon(i)
          lat_ori=list_lat(i)
        endif
      end do
      if( nx==-9999 )then
        print *, 'error: input area'
        stop
      endif

      allocate(nextx(nXX,nYY))
      allocate(catmx(nx,ny),catmy(nx,ny),flddif(nx,ny))

      allocate(flddph(nXX,nYY))
      allocate(flood(nx,ny))
! ===============================================
      fnextxy=trim(mapdir)//'nextxy.bin'
      fgridxy=trim(hires)//trim(area)//'.catmxy'
      fflddif=trim(hires)//trim(area)//'.flddif'

      fflood='./'//trim(area)//'.flood'

print *, fnextxy
      open(11, file=fnextxy, form='unformatted', access='direct', recl=2*nXX*nYY)
      read(11,rec=1) nextx
      close(11)

print *, fgridxy
      open(11, file=fgridxy, form='unformatted', access='direct', recl=2*nx*ny)
      read(11,rec=1) catmx
      read(11,rec=2) catmy
      close(11)

print *, fflddif
      open(11, file=fflddif, form='unformatted', access='direct', recl=4*nx*ny)
      read(11,rec=1) flddif
      close(11)

! =====

! modify to read flood depth from a file
      open(12, file=fflddph, form='unformatted', access='direct', recl=4*nXX*nYY)
      read(12,rec=trec) flddph
      close(12)

! sample, set flddph=5m
!      flddph(:,:)=-9999
!      do iYY=1, nYY
!        do iXX=1, nXX
!          if( nextx(iXX,iYY)/=-9999 )then
!            flddph(iXX,iYY)=5.0
!          endif
!        end do
!      end do
! ==========

      flood(:,:)=-9999
      do iy=1, ny
        do ix=1, nx
          if( catmx(ix,iy)>0 )then
            flood(ix,iy)=0
            iXX=catmx(ix,iy)
            iYY=catmy(ix,iy)
            if( flddph(iXX,iYY)>flddif(ix,iy) )then
              flood(ix,iy)=flddph(iXX,iYY)-flddif(ix,iy)
            endif
          endif
        end do
      end do

      open(11, file=fflood, form='unformatted', access='direct', recl=4*nx*ny)
      write(11,rec=1) flood
      close(11)

      end program proj_flood

