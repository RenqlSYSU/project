!***************************************************************************
!use the f_Qd, f_Qd_t, f_Qeddy and A to calc dzdt by SOR (Successive Over Relaxation)
!1. read forc and coe from binary file 
!2. calc dzdt by circulation 
!3. store the results
!
!                                            by Ql_Ren
!                                           2019/10/28
!******************************************************************************
program iteration 
    implicit none
    integer, parameter :: pr = 8 ,  nvar = 3, ntime = 37, nlev = 21 !ERA 
    integer, parameter :: nlat = 54 , nlon = 240 , ilog  = 10 , ifile = 12
    integer :: nc, nv, nt, nz, ny, nx, iter, irec 
    
    character(len=200) :: logname, fileout, filename
    character(len=2),dimension(nvar) :: var_name 
    
    real(kind=pr), parameter :: rf  = 1.25 ! relaxing factor
    real(kind=pr), parameter :: critical = 1e-8 
    real(kind=pr), parameter :: g  = 9.8 !M/(S*S)                     
    real(kind=pr), parameter :: cP = 1004.0 ! J/(K KG)  [ m2/(K s2) ] 
    real(kind=pr), parameter :: R  = 287.0                            
    real(kind=pr), parameter :: a  = 6378388 !the radius of earth, m  
    real(kind=pr), parameter :: pi = atan(1.0)*4                      

    real(kind=pr), dimension(nlon,nlat,nlev,ntime,nvar) :: forc, dzdt, dzdt0, term  
    real(kind=pr), dimension(nlon,nlat,nlev,ntime,nvar) :: coe111, coe110, coe101, coe121, coe011, coe211, f0 
    real(kind=pr), dimension(nlev) :: lev,dlev 
    real(kind=pr), dimension(nvar) :: diff

    lev  = (/1000,950,900,850,800,750,700,650,600,550,500,450,400,350,300,250,200,175,150,125,100/)
    dlev(2:nlev) = (- lev(1:(nlev-1)) + lev(2:nlev))*100 
    dlev(1) = dlev(2)
    
!    print*, "please input the number of case,1 or 2 or 3"
!    read(*,*) nc
     filename = "/home/ys17-19/renql/project/2020MiddleEastJet/data/ERA-Interim_DJF_3b6c.dat"
     fileout  = "/home/ys17-19/renql/project/2020MiddleEastJet/data/ERA-Interim_DJF_dzdt_re.dat"     
     logname  = "/home/ys17-19/renql/project/2020MiddleEastJet/f90_iter_infor_re-era.txt"
     var_name  = (/"turn","b1z ","b4  "/)  
     !fileout  = "/home/ys17-19/renql/project/2020MiddleEastJet/data/ERA-Interim_DJF_dzdt_advect_xy.dat"     
     !logname  = "/home/ys17-19/renql/project/2020MiddleEastJet/f90_iter_infor_advect-era.txt"
    ! var_name = (/"b1","b2","b3"/) !,"f_Qd_t "
    !var_name  = (/"b1x","b1y","b3x","b3y"/)  !this is the variable used to draw, can used as center string
    
    open(unit=ilog,file=logname,form='formatted',status='replace')
    write(ilog,*) "relaxing factor is ", rf
    write(ilog,*) "critical value is ", critical
    write(ilog,*) " "

    write(ilog,"(A200)") filename 

!==========================================================
!read the data
!============================================================
    open(unit=ifile,file=filename,form='binary',convert='LITTLE_ENDIAN',status='old',access='direct',recl=pr*nlat*nlon)
    irec = 1
    do nv = 1, nvar, 1
    do nt = 1, ntime, 1
    do nz = 1, nlev, 1
        read(ifile,rec=irec) ((forc(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
    end do
    end do
    end do
    
    do nv = 1, nvar, 1
    do nt = 1, ntime, 1
    do nz = 1, nlev, 1
        read(ifile,rec=irec) ((coe111(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((coe110(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((coe101(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((coe121(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((coe011(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((coe211(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
        read(ifile,rec=irec) ((f0(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
    end do
    end do
    end do
    close(ifile)
    write(ilog,*) "irec = ", irec 
    write(ilog,"(A20)") "Finish reading "
    write(ilog,"(A30,6(1x,E15.7))") "forc(5:10,5,5,1,1) are", forc(5:10,5,5,1,1) 
    write(ilog,"(A30,6(1x,E15.7))") "coe110(5:10,5,5) are", coe110(5:10,5,5,1,1) 

!==========================================================
!read the data
!============================================================
    dzdt = 0.0
    do iter = 1, 1000000, 1
        write(ilog,*) " "
        write(ilog,"(I5,A10)") iter, " iteration"
        
        dzdt0 = dzdt
        do nz = 2, nlev-1, 1
        do ny = 2, nlat-1, 1
        do nx = 2, nlon-1, 1
            term(nx,ny,nz,:,:) = coe111(nx,ny,nz,:,:)*dzdt(nx,ny,nz,:,:) + & 
                                 coe110(nx,ny,nz,:,:)*dzdt(nx+1,ny,nz,:,:) + coe121(nx,ny,nz,:,:)*dzdt(nx,ny+1,nz,:,:) + coe211(nx,ny,nz,:,:)*dzdt(nx,ny,nz+1,:,:) +&
                                 coe110(nx,ny,nz,:,:)*dzdt(nx-1,ny,nz,:,:) + coe101(nx,ny,nz,:,:)*dzdt(nx,ny-1,nz,:,:) + coe011(nx,ny,nz,:,:)*dzdt(nx,ny,nz-1,:,:) 
            dzdt(nx,ny,nz,:,:) = dzdt(nx,ny,nz,:,:) + (rf/coe111(nx,ny,nz,:,:))*(f0(nx,ny,nz,:,:)*a*a*forc(nx,ny,nz,:,:)-term(nx,ny,nz,:,:))
        end do
        end do
        end do
    
    !boundary conditions
        dzdt(1   ,:,:,:,:) = dzdt(nlon   ,:,:,:,:) !Periodic boundary for x
        !dzdt(1   ,:,:,:,:) = dzdt(2     ,:,:,:,:)
        !dzdt(nlon,:,:,:,:) = dzdt(nlon-1,:,:,:,:)
        
        dzdt(:,1   ,:,:,:) = 0 !dzdt(:,2     ,:,:,:)
        dzdt(:,nlat,:,:,:) = dzdt(:,nlat-1,:,:,:)
        
        dzdt(:,:,1   ,:,:) = 0  !lower boundary for A
        dzdt(:,:,nlev,:,:) = 0  !upper boundary for A
        
        do nv = 1, nvar ,1 
            write(ilog,*) "dzdt induced by "//trim(var_name(nv))//" is " ,dzdt(120,30,16,1,nv)
        end do
        do nv = 1, nvar ,1 
            diff(nv) = maxval(abs(dzdt(:,:,:,:,nv)-dzdt0(:,:,:,:,nv)))
            write(ilog,*) "Iterative difference of "//trim(var_name(nv))//" is " ,diff(nv)
        end do
        if (maxval(diff).lt.critical) then 
            write(ilog,*) "iteration has been finished"
            exit
        end if 
    end do
    write(ilog,*) dzdt(130:150:4,30,16,1,1)
    write(ilog,*) dzdt(130:150:4,30,10,1,1)
    
!========================================================================
!save the data
!=======================================================================
    open(unit=ifile,file=fileout,status='replace',form='binary',convert='LITTLE_ENDIAN',access='direct',recl=pr*nlat*nlon)
    irec = 1
    do nv = 1, nvar, 1
    do nt = 1, ntime, 1
    do nz = 1, nlev, 1
        write(ifile,rec=irec) ((dzdt(nx,ny,nz,nt,nv),nx=1,nlon),ny=1,nlat)
        irec = irec + 1
    end do
    end do
    end do
    write(ilog,*) "irec = ", irec
    close(ifile)
    close(ilog)
end program

