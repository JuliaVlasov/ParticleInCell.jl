!UA scheme for 4d Vlasov in Fluid-scaling with b(x)
!External E
!finite difference solver
!with new IMEX1 and IMEX2
program efd

use fft_m
use mesh_fields_m, only: fields_2d_t, mesh_t, init_mesh, init_fields
use particles_m, only: init_particles_2d, particles_t
use m6_compute_rho_m, only: compute_rho_m6_real

implicit none

type(mesh_t)      :: mesh
type(fields_2d_t) :: f
type(particles_t) :: p

integer, parameter :: ntau = 16
integer, parameter :: npp  = 204800

real(8)    :: xxt(2)
real(8)    :: bx
real(8)    :: ds
real(8)    :: interv
real(8)    :: dtau
real(8)    :: tau(ntau)
real(8)    :: ltau(ntau)
real(8)    :: energy
real(8)    :: Et(2,ntau)
real(8)    :: ave(2)
real(8)    :: ave2(2)
real(8)    :: ave3(2)

complex(8) :: pl(       ntau)
complex(8) :: ql(       ntau)
complex(8) :: xf( ntau,2)
complex(8) :: yf( ntau,2)
complex(8) :: gx(ntau,2)
complex(8) :: gy(ntau,2)
complex(8) :: temp(   ntau,2)
complex(8) :: tilde(  ntau,2)
complex(8) :: h(      ntau,2)
complex(8) :: r(      ntau,2)
complex(8) :: fx(     ntau,2)
complex(8) :: fy(     ntau,2)
complex(8) :: xt(     ntau,2)
complex(8) :: yt(     ntau,2)

real(8) :: time
real(8) :: xmin
real(8) :: xmax
real(8) :: ymin
real(8) :: ymax
integer :: istep
integer :: iplot
integer :: iargc
integer :: n,m
integer :: i
integer :: j
integer :: error

real(8) :: aux1, aux2

character(len=272) :: argv

integer, parameter :: nx = 128
integer, parameter :: ny = 64

integer(8) :: nbpart
integer    :: nstep

real(8) :: alpha
real(8) :: dt
real(8) :: dx
real(8) :: dy
real(8) :: ep
real(8) :: kx
real(8) :: ky
real(8) :: dimx
real(8) :: dimy
real(8) :: tfinal
real(8) :: pi
real(8) :: poids
real(8) :: x1, x2, v1, v2

complex(8), parameter :: im = (0d0, 1d0)


pi    = 4d0 * atan(1d0)
alpha = 0.05d0 
kx    = 0.50_8
ky    = 1.0d0   
dimx  = 2*pi/kx
dimy  = 2*pi/ky  
poids = dimx * dimy 

dx = dimx / nx
dy = dimy / ny

dt     = pi/2.d0/(2.0d0**3)!Used to determine N in MRC, not the macro-step
tfinal = pi/2.d0
nstep  = nint(tfinal/dt)

time  = 0.d0

ep   = 0.001d0
dtau = 2.0d0*pi/ntau

call init_fft( ntau, temp(:,1), tilde(:,1) )

m = ntau/2
ltau=(/ (n, n=0,m-1), (n, n=-m,-1 )/)
do i=1,ntau
  tau(i) = (i-1)*dtau
end do

xmin = 0.0_8; xmax = dimx
ymin = 0.0_8; ymax = dimy

nbpart = npp

call init_mesh( mesh, xmin, xmax, nx, ymin, ymax, ny )
call init_fields( f, mesh )
call init_particles_2d( p, nbpart, mesh, alpha, kx )

print"('ep = ', g15.3)", ep
print"('dt = ', g15.3)", dt

print*, sum(p%x(1,:)), sum(p%x(2,:)), sum(p%v(1,:)), sum(p%v(2,:))

do m=1,1!nbpart
    
    x1 = p%x(1,m)
    x2 = p%x(2,m)
    v1 = p%v(1,m)
    v2 = p%v(2,m)

    time  = 0.d0
    bx    = 1.d0+0.5d0*sin(x1)*sin(x2)
    ds    = dt*bx

    pl(1) = ds
    ql(1) = ds**2/2.0d0

    do i=2,ntau
        pl(i) = ep*im*(exp(-im*ltau(i)*ds/ep)-1.0d0)/ltau(i)
        ql(i) = ep*(ep*(1.0d0-exp(-im*ltau(i)*ds/ep))-im*ltau(i)*ds)/ltau(i)**2
    end do

    xt(:,1) = x1
    xt(:,2) = x2

    yt(:,1) = v1
    yt(:,2) = v2

    !--preparation initial data--
    temp(1,1) = v1/bx
    temp(1,2) = v2/bx
    do n=1,ntau
        h(n,1)=ep*(sin(tau(n))*temp(1,1)-cos(tau(n))*temp(1,2))
        h(n,2)=ep*(sin(tau(n))*temp(1,2)+cos(tau(n))*temp(1,1))
    end do

    xt(:,1)=x1+h(:,1)-h(1,1)
    xt(:,2)=x2+h(:,2)-h(1,2)!xt1st

    p%e(1,m)=(0.5d0*cos(x1/2.d0)*sin(x2))*(1.d0+0.5d0*sin(time))
    p%e(2,m)=(sin(x1/2.d0)*cos(x2))*(1.d0+0.5d0*sin(time))

    do n=1,ntau
        interv=(1.d0+0.5d0*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx
        r(n,1) =  interv*v2
        r(n,2) = -interv*v1
    end do
    print*, r(1,1), r(1,2)

    print*, sum(r)
    call fft(r, tilde)
    print*, sum(tilde)
    ave = real(tilde(1,:))/ep!dot{Y}_underline^0th

    print*, " ave : ", ave
    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do

    tilde(1,:)=0.0d0
    call ifft( tilde, r)

    do n=1,ntau
        r(n,1)=ep*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))/bx+r(n,1)
        r(n,2)=ep*(sin(tau(n))*p%e(2,m)-cos(tau(n))*p%e(1,m))/bx+r(n,2)
    end do

    yt(:,1)=v1+(r(:,1)-r(1,1))
    yt(:,2)=v2+(r(:,2)-r(1,2))!yt1st

    print*, " yt : ",  sum(yt(:,1)), sum(yt(:,2))
    !--more preparation--

    do n=1,ntau
        temp(n,1)=ep*(cos(tau(n))*yt(n,1)+sin(tau(n))*yt(n,2))/bx
        temp(n,2)=ep*(cos(tau(n))*yt(n,2)-sin(tau(n))*yt(n,1))/bx
    end do

    call fft(temp, tilde)

    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do

    tilde(1,:)=0.0d0
    call ifft(tilde, h)

    do n=1,ntau
        h(n,1)=h(n,1)-ep**2/bx*(-cos(tau(n))*ave(1)-sin(tau(n))*ave(2))
        h(n,2)=h(n,2)-ep**2/bx*(-cos(tau(n))*ave(2)+sin(tau(n))*ave(1))!h2nd
    end do

    xt(:,1)=x1+h(:,1)-h(1,1)
    xt(:,2)=x2+h(:,2)-h(1,2)!x2nd,residue O(eps^3)

    print*, " xt : ",  sum(xt(:,1)), sum(xt(:,2))

    p%e(1,m)=(0.5d0*cos(x1/2.d0)*sin(x2))*0.5d0*cos(time)
    p%e(2,m)=(sin(x1/2.d0)*cos(x2))*0.5d0*cos(time)!partial_tE

    print*, p%e(1,1), p%e(2,1)

    do n=1,ntau

        interv=(1.d0+0.5d0*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx

        fx(n,1)=interv*ave(2)
        fx(n,2)=-interv*ave(1)

        fy(n,1)=ep/bx*(sin(tau(n))*ave(1)-cos(tau(n))*ave(2))
        fy(n,2)=ep/bx*(cos(tau(n))*ave(1)+sin(tau(n))*ave(2))!partial_sh^1st

        interv=cos(x1)*sin(x2)*fy(n,1)+sin(x1)*cos(x2)*fy(n,2)

        fy(n,1)=interv/bx/2.d0*v2+fx(n,1)
        fy(n,2)=-interv/bx/2.d0*v1+fx(n,2)

        fx(n,1)=ep/bx**2*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        fx(n,2)=ep/bx**2*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))

    end do

    temp = fy + fx
    print*, " fx   : ", sum(fx(:,1)), sum(fx(:,2))
    print*, " fy   : ", sum(fy(:,1)), sum(fy(:,2))
    print*, " temp : ", sum(temp)
    print*, ' ep   : ', ep
    print*, ' bx   : ', bx
    print*, ' ave  : ', ave
    stop

    call fft(temp, tilde)

    do n=2,ntau
        fx(n,:)=-im*tilde(n,:)/ltau(n)
        tilde(n,:)=-tilde(n,:)/ltau(n)**2
    end do

    fx(1,:)    = 0.0d0
    tilde(1,:) = 0.0d0

    call ifft(tilde, temp)

    r = - ep * temp

    call ifft(fx, fy)

    gx=fy!partial_sr^1st!!!
    do n=1,ntau
        tilde(n,1)=(cos(tau(n))*fy(n,1)+sin(tau(n))*fy(n,2))/bx
        tilde(n,2)=(cos(tau(n))*fy(n,2)-sin(tau(n))*fy(n,1))/bx
    end do

    call fft(tilde, temp)

    ave2=temp(1,:)/ntau!ddot{X_underline}!!!!
    do n=1,ntau

        p%e(1,m)=(0.5d0*cos(real(xt(n,1))/2.d0)*sin(real(xt(n,2))))*(1.d0+0.5d0*sin(time))
        p%e(2,m)=(sin(real(xt(n,1))/2.d0)*cos(real(xt(n,2))))*(1.d0+0.5d0*sin(time))

        interv=(1.d0+0.5d0*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx

        temp(n,1)=interv*yt(n,2)+ep/bx*(-sin(tau(n))*p%e(2,m)+cos(tau(n))*p%e(1,m))
        temp(n,2)=-interv*yt(n,1)+ep/bx*(sin(tau(n))*p%e(1,m)+cos(tau(n))*p%e(2,m))
    end do

    call fft( temp, tilde)

    xf(1,:)=tilde(1,:)/ep!dot{Y_underline}
    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do
    tilde(1,:)=0.0d0
    call ifft(tilde, temp)
    r=r+temp

    yt(:,1)=v1+r(:,1)-r(1,1)
    yt(:,2)=v2+r(:,2)-r(1,2)

    !----end more preparation
    !--even more preparation--

    do n=1,ntau
        temp(n,1)=(cos(tau(n))*r(n,1)+sin(tau(n))*r(n,2))/bx
        temp(n,2)=(cos(tau(n))*r(n,2)-sin(tau(n))*r(n,1))/bx
    end do

    call fft( temp, tilde)

    temp(1,:)=tilde(1,:)
    ave3=temp(1,:)!dot{X_underline}!!!
    interv=cos(x1)*sin(x2)*temp(1,1)+sin(x1)*cos(x2)*temp(1,2)

    fx(1,1)=interv/ep/bx*v2/2.d0
    fx(1,2)=-interv/ep/bx*v1/2.d0

    do n=1,ntau
        temp(n,1)=(1.d0+0.5d0*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx
    end do

    call fft(temp(:,1), tilde(:,1))

    fx(1,1)=fx(1,1)+tilde(1,1)/ep*ave(2)
    fx(1,2)=fx(1,2)-tilde(1,1)/ep*ave(1)!ddot{Y}_underline^0th!!!!

    do n=1,ntau
        yf(n,:)=xf(1,:)+fy(n,:)
        temp(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))
        temp(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))
    end do

    call fft( temp, tilde)

    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do

    tilde(1,:)=0.0d0

    call ifft( tilde, temp)

    do n=1,ntau
        fy(n,1)=temp(n,1)*ep/bx-ep**2/bx*(-cos(tau(n))*fx(1,1)-sin(tau(n))*fx(1,2))
        fy(n,2)=temp(n,2)*ep/bx-ep**2/bx*(-cos(tau(n))*fx(1,2)+sin(tau(n))*fx(1,1))!partial_sh2!!!
    end do

    call fft( fy, tilde)

    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do
    tilde(1,:)=0.0d0

    call ifft( tilde, temp)

    h=-ep*temp

    do n=1,ntau
        temp(n,1)=ep*(cos(tau(n))*yt(n,1)+sin(tau(n))*yt(n,2))/bx
        temp(n,2)=ep*(cos(tau(n))*yt(n,2)-sin(tau(n))*yt(n,1))/bx
    end do

    call fft( temp, tilde)

    do n=2,ntau
        tilde(n,:)=-im*tilde(n,:)/ltau(n)
    end do

    tilde(1,:)=0.0d0

    call ifft( tilde, temp)

    h = h + temp

    xt(:,1)=x1+h(:,1)-h(1,1)
    xt(:,2)=x2+h(:,2)-h(1,2)
    
    !--end even more
    !--iteration
    do istep = 1, nstep

        !---imex2 New---
        call compute_fy( xt, yt)

        gy=yt+ds/2.d0*fy

        call fft( gy, fy)

        do n=1,ntau
            fy(n,:)=fy(n,:)/(1.0d0+im*ds/2.d0*ltau(n)/ep)
        end do

        call ifft( fy, yf)!yt(tn+1/2)

        do n=1,ntau
            fx(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))/bx
            fx(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))/bx
        end do

        gx = xt + ds/2.d0 * fx

        call fft( gx, fx)

        do n=1,ntau
            fx(n,:)=fx(n,:)/(1.0d0+im*ds/2.d0*ltau(n)/ep)
        end do

        call ifft( fx, xf)!xt(tn+1/2)

        time = time + dt/2.d0

        call compute_fy( xf, yf)

        call fft( fy, gy)
        call fft( yt, yf)

        do n=1,ntau
            fy(n,:)=(yf(n,:)*(1.0d0-im*ds/ep/2.0d0*ltau(n)) &
&            +ds*gy(n,:))/(1.0d0+im*ds/2.0d0*ltau(n)/ep)
        end do

        yf=yt

        call ifft( fy, yt)!yt(tn+1)

        yf=(yt+yf)/2.d0
        do n=1,ntau
            fx(n,1)=(cos(tau(n))*yf(n,1)+sin(tau(n))*yf(n,2))/bx
            fx(n,2)=(cos(tau(n))*yf(n,2)-sin(tau(n))*yf(n,1))/bx
        end do

        call fft( fx, gx)
        call fft( xt, xf)

        do n=1,ntau
            fx(n,:)=(xf(n,:)*(1.0d0-im*ds/ep/2.0d0*ltau(n)) &
        &       +ds*gx(n,:))/(1.0d0+im*ds/2.0d0*ltau(n)/ep)
        end do

        call ifft( fx, xt)!xt(tn+1)

        time=time+dt/2.d0

        !---end imex2 New---

    end do

    call fft( xt,tilde)

    temp(1,:)=0.d0
    do n=1,ntau
        temp(1,:)=temp(1,:)+tilde(n,:)*exp(im*ltau(n)*tfinal*bx/ep)
    end do

    xxt=real(temp(1,:))

    call apply_bc()

    p%x(1,m) = xxt(1)
    p%x(2,m) = xxt(2)

    call fft(yt,tilde)

    temp(1,:)=0.d0
    do n=1,ntau
        temp(1,:)=temp(1,:)+tilde(n,:)*exp(im*ltau(n)*tfinal*bx/ep)
    end do

    p%v(1,m)=real(cos(tfinal*bx/ep)*temp(1,1)+sin(tfinal*bx/ep)*temp(1,2))
    p%v(2,m)=real(cos(tfinal*bx/ep)*temp(1,2)-sin(tfinal*bx/ep)*temp(1,1))

end do
print*, sum(p%v(1,:))+857.95049281063064, sum(p%v(2,:))+593.40700170710875 


call compute_rho_m6_real(f, p)
stop
open(unit=851,file='UA3ep0001T1.dat')
do i=1,nx
do j=1,ny
write(851,*)i*dx, j*dy, f%rho(i,j)
end do
write(851,*)
end do
close(851)
!call calcul_energy( p, f,energy )
!open(unit=851,file='UA3ep0001T1v.dat')
!do i=1,nx
!do j=1,ny
!write(851,*)i*dx, j*dy, f%rho(i,j)
!end do
!write(851,*)
!end do
!close(851)

contains


subroutine compute_fy( xt, yt )

    complex(8) :: xt(:,:)
    complex(8) :: yt(:,:)

    do n=1,ntau
        Et(1,n)=(0.5d0*cos(real(xt(n,1))/2.d0)*sin(real(xt(n,2))))*(1.d0+0.5d0*sin(time))
        Et(2,n)=(cos(real(xt(n,2)))*sin(real(xt(n,1))/2.d0))*(1.d0+0.5d0*sin(time))
        interv=(1.d0+0.5d0*sin(real(xt(n,1)))*sin(real(xt(n,2)))-bx)/bx/ep
        temp(n,1)=(cos(tau(n))*Et(1,n)-sin(tau(n))*Et(2,n))/bx
        temp(n,2)=(cos(tau(n))*Et(2,n)+sin(tau(n))*Et(1,n))/bx
        fy(n,1)=temp(n,1)+interv*yt(n,2)
        fy(n,2)=temp(n,2)-interv*yt(n,1)
    end do

end subroutine compute_fy

subroutine apply_bc()

    do while ( xxt(1) > xmax )
        xxt(1) = xxt(1) - dimx
    end do

    do while ( xxt(1) < xmin )
        xxt(1)= xxt(1) + dimx
    end do

    do while ( xxt(2) > ymax )
        xxt(2)  = xxt(2)  - dimy
    end do

    do while ( xxt(2)  < ymin )
        xxt(2) = xxt(2)  + dimy
    end do

end subroutine apply_bc


end program efd



