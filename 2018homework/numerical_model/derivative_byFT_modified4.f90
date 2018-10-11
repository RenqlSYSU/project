
PROGRAM SPECTRAL
  IMPLICIT NONE
  
  REAL(KIND=4), PARAMETER :: PI  = 4.*atan(1.), &
                             PI2 = 2.*PI

  REAL(KIND=4), DIMENSION(:,:),   ALLOCATABLE :: x,f,ft,df
  COMPLEX(KIND=4), DIMENSION(:,:),ALLOCATABLE :: t
  REAL, DIMENSION(3) :: err_rms, err_rms_ref,err_inf
  INTEGER(KIND=4) n,mode
  REAL(KIND=4) x0,x1
  PARAMETER(x0=-1., x1=1.)
  INTEGER ni,i,j

  WRITE(*,*) 'Usage: ./Derivative_via_FT.x <n> <mode>'
  WRITE(*,*) '       <n>    Positive Integer Number of collocation points'
  WRITE(*,*) '       <mode> Function to test with'
  WRITE(*,*) '              1: y = sin(pi*x) '
  WRITE(*,*) '              2: y = exp(-5x^2)'
  WRITE(*,*) '              3: y = 1-abs(x)  ' 
  WRITE(*,*) 'Using Default values: n=128    '
  n=128
  ALLOCATE(x(3,n),f(3,n),ft(3,n),t(3,n),df(3,n))

  CALL INIT(n,x0,x1,x,f,df,mode)
  
  do j=1,3
  CALL FFT_DERIVATIVE(n,x0,x1,f(j,:),t(j,:),ft(j,:))
  enddo

  ni = 30
  WRITE(*,*) '                x: ',x(0,ni)
  WRITE(*,*) '1: y = sin(pi*x) = ',f(1,ni),df(1,ni),ft(1,ni)
  WRITE(*,*) '2: y = exp(-5x^2)= ',f(2,ni),df(2,ni),ft(2,ni)
  WRITE(*,*) '3: y = 1-abs(x)  = ',f(3,ni),df(3,ni),ft(3,ni) 
  ni = 100
  WRITE(*,*) '                x: ',x(0,ni)
  WRITE(*,*) '1: y = sin(pi*x) = ',f(1,ni),df(1,ni),ft(1,ni)
  WRITE(*,*) '2: y = exp(-5x^2)= ',f(2,ni),df(2,ni),ft(2,ni)
  WRITE(*,*) '3: y = 1-abs(x)  = ',f(3,ni),df(3,ni),ft(3,ni) 
  
  err_rms =0.
  err_rms_ref=0.
  do j=1,3
  err_inf(j) = MAXVAL(ABS(ft(j,:)-df(j,:)))
  end do

  DO i=1,n
     err_rms = err_rms + (ft(:,i)-df(:,i))**2
     err_rms_ref = err_rms_ref+ df(:,i)**2
  ENDDO

  WRITE(*,*) 'INFINITY NORM ERROR: ',err_inf
  WRITE(*,*) 'RMS ERROR (ABSOLUTE):',err_rms/n
  WRITE(*,*) 'RMS ERROR (RELATIVE):',err_rms/err_rms_ref

  OPEN(20,FILE='derivative4.dat',FORM='binary',ACCESS='DIRECT',RECL=4*3*N)
  WRITE(20,rec=1) ((x(i,j),j=1,n),i=1,3) 
  WRITE(20,rec=2) ((f(i,j),j=1,n),i=1,3)
  WRITE(20,rec=3) ((df(i,j),j=1,n),i=1,3)  !real first derivative
  WRITE(20,rec=4) ((ft(i,j),j=1,n),i=1,3)  !first derivative by fourier transform
  CLOSE(20)

  STOP 'SUCCESS'

CONTAINS
  SUBROUTINE INIT(n,x0,x1,x,f,df,mode)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN)   ::  n,mode
    REAL(KIND=4),    INTENT(IN)   ::  x0,x1
    REAL(KIND=4),    DIMENSION(3,n) ::  x,f,df
    INTEGER i
    REAL dx

    dx = (x1-x0)/n
    DO i=1,n
       x(:,i) = x0 + (i-1)*dx

       f(1,i) = SIN(1*PI*x(1,i))
       df(1,i)= PI*COS(1*PI*x(1,i))

       f(2,i) = exp(-5*x(2,i)**2)
       df(2,i)= -10*x(2,i)*exp(-5*x(2,i)**2)

       f(3,i) = 1-abs(x(3,i))
    if(x(3,i).eq.0) then
       df(3,i)= -1
    else
       df(3,i)= -abs(x(3,i))/x(3,i)
    endif
    ENDDO
  END SUBROUTINE INIT

  SUBROUTINE FFT_DERIVATIVE(n,x0,x1,f,t,df)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    REAL(KIND=4),    DIMENSION(n), INTENT(IN)   :: f
    REAL(KIND=4),    DIMENSION(n), INTENT(OUT)  :: df
    REAL(KIND=4),    INTENT(IN) :: x0,x1
    COMPLEX(KIND=4), DIMENSION(n), INTENT(INOUT):: t
    REAL k

    CALL FFT_F(n,f,t)
    DO i=1,n
       IF ( i-1<=n/2) THEN
          k=i-1
       ELSE
          k=-(n-i+1)
       ENDIF
       t(i) = t(i)*CMPLX(0,k)*(2*PI/(x1-x0))
    ENDDO
    CALL FFT_B(n,t,df)
  END SUBROUTINE FFT_DERIVATIVE

  SUBROUTINE FFT_F(n,r,c)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    REAL(KIND=4),    DIMENSION(n), INTENT(IN) :: r
    COMPLEX(KIND=4), DIMENSION(n), INTENT(OUT):: c
    INTEGER j,k

    DO k=1,n
       c(k) = CMPLX(0.,0.)
       DO j=1,n
          c(k) = c(k) + r(j) * EXP( CMPLX(0, -PI2*(j-1)*(k-1)/n) )
       ENDDO
       c(k) = c(k) / SQRT(REAL(n))
    ENDDO
  END SUBROUTINE FFT_F


  SUBROUTINE FFT_B(n,c,r)
    IMPLICIT NONE
    INTEGER(KIND=4), INTENT(IN) :: n
    REAL(KIND=4),    DIMENSION(n), INTENT(OUT) :: r
    COMPLEX(KIND=4), DIMENSION(n), INTENT(IN):: c
    COMPLEX(KIND=4) :: c_loc
    INTEGER j,k,k_loc

    DO j=1,n
       c_loc = 0.
       DO k=1,n
          c_loc = c_loc + c(k) * EXP( CMPLX(0, PI2*(j-1)*(k-1)/n) )
       ENDDO
       r(j) = REAL(c_loc) / SQRT(REAL(n))
    ENDDO
  END SUBROUTINE FFT_B
END PROGRAM SPECTRAL
