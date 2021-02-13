module particules

use iso_c_binding

use zone
use quietstart

implicit none

integer, private :: ipart 
integer, private :: i, j

CONTAINS

subroutine plasma( p )

real(c_double), allocatable :: p(:,:)
real(c_double) :: speed, theta, vth, n
real(c_double) :: a, b, eps, R
integer :: k

eps = 1.d-12

vth =  1.
nbpart = 100*(nx)*(ny)
n = 1.d0/nbpart

allocate(p(7, nbpart))

do k=0,nbpart-1

   speed = vth * sqrt(-2 * log( (k+0.5)*n ))

   theta = trinary_reversing( k ) * 2 * pi

   a = 0; b = dimx ! 2*pi/kx 
   R = bit_reversing( k )
   call dichotomie_x(a,b,R,eps) 
   p(1,k+1) = a
   p(2,k+1) = dimy * penta_reversing( k ) 
   p(3,k+1) = speed * cos(theta)
   p(4,k+1) = speed * sin(theta)

enddo

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpolation( f, p ) bind(C, name="interpolation")

real(c_double) :: p(:,:)
real(c_double) :: f(:,:,:)
real(c_double) :: a1, a2, a3, a4
real(c_double) :: xp, yp, dxp, dyp
integer :: ip1, jp1

do ipart=1,nbpart

   xp = p(1,ipart) / dx
   yp = p(2,ipart) / dy

   i = floor( xp ) + 1
   j = floor( yp ) + 1

   dxp = xp - i + 1
   dyp = yp - j + 1

   a1 = (1-dxp) * (1-dyp) 
   a2 = dxp * (1-dyp) 
   a3 = dxp * dyp 
   a4 = (1-dxp) * dyp 

   ip1 = mod1(i+1,nx)
   jp1 = mod1(j+1,ny)

   p(5:7,ipart) = a1 * f(1:3,i,j) + a2 * f(1:3,ip1,j) &
                + a3 * f(1:3,ip1,jp1) + a4 * f(1:3,i,jp1) 

end do

end subroutine interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine push_v( p ) bind(C, name="push_v")

real(c_double) :: p(:,:)
real(c_double) :: dum
real(c_double) :: tantheta, sintheta
real(c_double) :: hdt, v1, v2, e1, e2, b3

hdt = 0.5 * dt

do ipart = 1, nbpart

   v1 = p(3,ipart)
   v2 = p(4,ipart)
   e1 = p(5,ipart)
   e2 = p(6,ipart)
   b3 = p(7,ipart)

   v1 = v1 + hdt * e1
   v2 = v2 + hdt * e2

   tantheta = hdt * b3
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   v1 = v1 + v2 * tantheta
   v2 = v2 - v1 * sintheta
   v1 = v1 + v2 * tantheta

   p(3,ipart) = v1 + hdt * e1
   p(4,ipart) = v2 + hdt * e2

end do

end subroutine push_v

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine push_x( p, coef ) bind(C, name="push_x")

real(c_double) :: p(:,:)
real(c_double) :: coef

do ipart=1,nbpart
   p(1,ipart) = p(1,ipart) + p(3,ipart) * dt * coef
   p(2,ipart) = p(2,ipart) + p(4,ipart) * dt * coef
   p(1,ipart) = modulo(p(1,ipart), dimx)
   p(2,ipart) = modulo(p(2,ipart), dimy)
end do   

end subroutine push_x

subroutine deposition( p, f0, f1 ) bind(C, name="deposition")

real(c_double) :: p(:,:), f0(:,:,:), f1(:,:,:)
real(c_double) :: a1, a2, a3, a4, dum, xp, yp, dxp, dyp
integer :: ip1, jp1

f1(4:5,:,:) = 0.d0

do ipart=1,nbpart

   xp = p(1,ipart) / dx
   yp = p(2,ipart) / dy

   i = floor( xp ) + 1
   j = floor( yp ) + 1

   ip1 = mod1(i+1,nx)
   jp1 = mod1(j+1,ny)

   dxp = xp - i + 1
   dyp = yp - j + 1

   a1 = (1-dxp) * (1-dyp) 
   a2 = dxp * (1-dyp) 
   a3 = dxp * dyp 
   a4 = (1-dxp) * dyp 

   dum = p(3,ipart) / (dx*dy) * dimx * dimy / nbpart

   f1(1,i,j)     = f1(1,i,j)     + a1*dum  
   f1(1,ip1,j)   = f1(1,ip1,j)   + a2*dum 
   f1(1,ip1,jp1) = f1(1,ip1,jp1) + a3*dum 
   f1(1,i,jp1)   = f1(1,i,jp1)   + a4*dum 

   dum = p(4,ipart) / (dx*dy) * dimx * dimy / nbpart

   f1(2,i,j)     = f1(2,i,j)     + a1*dum  
   f1(2,ip1,j)   = f1(2,ip1,j)   + a2*dum 
   f1(2,ip1,jp1) = f1(2,ip1,jp1) + a3*dum 
   f1(2,i,jp1)   = f1(2,i,jp1)   + a4*dum 

end do

do i=1,nx
do j=1,ny+1
   f0(1,i,j) = 0.5 * (f1(1,i,mod1(j,ny))+f1(1,mod1(i+1,nx),mod1(j,ny)))
end do
end do

do i=1,nx+1
do j=1,ny
   f0(2,i,j) = 0.5 * (f1(2,mod1(i,nx),j)+f1(2,mod1(i,nx),mod1(j+1,ny)))
end do
end do

end subroutine deposition

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
