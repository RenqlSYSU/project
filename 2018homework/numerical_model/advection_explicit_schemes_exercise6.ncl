;************************************************
;advection_explicit_schemes_exercise6.ncl
;
;This is the sixth homework of Atmospheric Modeling
;FTCS, Lax, and Lax-Wendroff schemes are used to difference the 
;advection equation. Then plots the result.
;
;                                       by Renql
;                                      2018/11/29
;************************************************

begin
;---Define parameters------------------
c  = 1.0    ;phase velocity
dx = 0.02   ;grid spacing
dt = 0.002  ;time spacing
lx = 1.0    ;the length of x
nx = 50     ;the number of grids
nstp = lx/(c*dt)+1    ;integration steps

;---Boundary and initial value--------------
k     = 1.0
theta = 1.0
x0    = 0
f = new((/nstp,nx/),float)
x = fspan(-(lx-dx)/2, (lx-dx)/2, nx)
f(0,:) = cos(k*(x-x0))*exp(-(x-x0)^2/(2*theta^2))
print(x(0:10)) 
print(f(0,0:10))

;---Explicit schemes for the advection equation---------
opt = 0
if(opt.eq.0) then
    title = "FTCS"
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f1(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,ix+1) - f(istp-1,nx-1))
        else if(ix.eq.(nx-1)) then
            f1(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,0) - f(istp-1,ix-1))
        else
            f1(istp,ix) = f(istp-1,ix) - c*dt/2/dx*( f(istp-1,ix+1) - f(istp-1,ix-1))
        end if
    end do
    end do
else if(opt.eq.1) then
    title = "Lax"
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f1(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,ix+1) + (0.5 + c*dt/2/dx)*f(istp-1,nx-1)
        else if(ix.eq.(nx-1)) then
            f1(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,0) + (0.5 + c*dt/2/dx)*f(istp-1,ix-1)
        else
            f1(istp,ix) = (0.5 - c*dt/2/dx)*f(istp-1,ix+1) + (0.5 + c*dt/2/dx)*f(istp-1,ix-1)
        end if
    end do
    end do
else
    title = "Lax-Wendroff"
    do istp = 1, nstp-1, 1
    do ix = 0, nx-1, 1
        if(ix.eq.0) then
            f1(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,ix+1) + (c*dt/dx+1)*f(istp-1,nx-1))
        else if(ix.eq.(nx-1)) then
            f1(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,0) + (c*dt/dx+1)*f(istp-1,ix-1))
        else
            f1(istp,ix) = (1-(c*dt/dx)^2)*f(istp-1,ix) + c*dt/2/dx*( (c*dt/dx-1)*f(istp-1,ix+1) + (c*dt/dx+1)*f(istp-1,ix-1))
        end if
    end do
    end do
end if


