program read_input_jet
   implicit none

   integer, parameter :: nz = 64
   integer, parameter :: ny = 80 
   real, dimension(nz,ny) :: u,r,t,zk
   integer :: ny_in, nz_in, j,k
   real, dimension(ny,nz) :: field_in

! this code assumes it is called on processor 0 only

   OPEN(unit=10, file='input_jet', form='unformatted', status='old')
   REWIND(10) !Back to the head of the file
   read(10) ny_in,nz_in
   write(*,*) ny_in, nz_in

   if((ny_in /= ny ) .or. (nz_in /= nz)) then
     write(*,*) ' error in input jet dimensions '
     write(*,*) ' ny, ny_input, nz, nz_input ', ny, ny_in, nz,nz_in
     write(*,*) ' error exit '
     write(*,*) ' error in input jet dimensions ' 
   end if
   
   read(10) field_in
   do j=1,ny
   do k=1,nz
     u(k,j) = field_in(j,k)
   enddo
   enddo
  
   read(10) field_in
   do j=1,ny
   do k=1,nz
     t(k,j) = field_in(j,k)
   enddo
   enddo

   read(10) field_in
   do j=1,ny
   do k=1,nz
     r(k,j) = field_in(j,k)
   enddo
   enddo

   do j=1,ny
   do k=1,nz
     zk(k,j) = 125. + 250.*float(k-1)
   enddo
   enddo
   close(10)

   OPEN(unit=11, file='input_jet.txt')
   write(11,*) 'u'
   do k=1,nz
   write(11,*) u(k,:)
   enddo

   write(11,*) 'temperature'
   do k=1,nz
   write(11,*) t(k,:)
   enddo

   write(11,*) 'density'
   do k=1,nz
   write(11,*) r(k,:)
   enddo
   
   write(11,*) 'height'
   do k=1,nz
   write(11,*) zk(k,:)
   enddo

   close(11)

end

