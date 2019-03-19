      program vec2map
! ================================================
      implicit none
! === Basic Parameters ===========================
! Resolution
      integer             ::  nx                  !! number of grids in horizontal
      integer             ::  ny                  !! number of grids in vertical
      integer             ::  nseqmax
      integer             ::  nflp
! Input
      character*256       ::  outdir, varname, cyear, buf
      integer             ::  iyear, regionall
! Index
      integer             ::  ix, iy, iseq
      integer             ::  irec, nrec
      integer             ::  iregion
! Var
      integer,allocatable ::  seqx(:),seqy(:)
      real,allocatable    ::  varmap(:,:), varvec(:)
!
      real                ::  rmis
      parameter              (rmis=1.e20)
! File
      character*256       ::  fmap, fvec, fseq, fmapdata
      character*256       ::  csufbin, csufvec
      parameter              (csufbin='.bin')
      parameter              (csufvec='.vec')
      character*256       ::  seqname
      parameter              (seqname='ind_xy')
! Function
      integer             ::  ndays
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
      read(11,*) buf, regionall
      read(11,*) buf, nseqmax
      close(11)

      allocate(varmap(nx,ny))

      allocate(seqx(nseqmax))
      allocate(seqy(nseqmax))
      allocate(varvec(nseqmax))

      nrec=ndays(iyear)

      fmap=trim(outdir)//trim(varname)//trim(cyear)//trim(csufbin)
      open(21,file=fmap,form='unformatted',access='direct',recl=4*nx*ny)

      do irec=1, nrec

        varmap(:,:)=rmis

        do iregion=1, regionall
          if( regionall>=2 )then
            write(buf,'(a1,i2.2)') '-', iregion
            fvec=trim(outdir)//trim(varname)//trim(cyear)//trim(buf)//trim(csufvec)
            fseq=trim(outdir)//trim(seqname)//trim(buf)//trim(csufvec)
          else
            fvec=trim(outdir)//trim(varname)//trim(cyear)//trim(csufvec)
            fseq=trim(outdir)//trim(seqname)//trim(csufvec)
          endif

          open(11,file=fseq,form='unformatted',access='direct',recl=4*nseqmax)
          read(11,rec=1) seqx
          read(11,rec=2) seqy
          close(11)

          open(11,file=fvec,form='unformatted',access='direct',recl=4*nseqmax,status='old')
          read(11,rec=irec) varvec
          close(11)

          do iseq=1, nseqmax
            if( seqx(iseq)>0 )then
              ix=seqx(iseq)
              iy=seqy(iseq)
              varmap(ix,iy)=varvec(iseq)
            endif
          end do
        end do

        write(21,rec=irec) varmap

      end do

      close(21)

      end program vec2map



!!==================================================
      integer function ndays(iyear)
      implicit none
      integer :: iyear

      ndays=365
      if ( mod(iyear,400)==0 .or. (mod(iyear,100)/=0 .and. mod(iyear,4)==0 )) ndays=366

      end function ndays


