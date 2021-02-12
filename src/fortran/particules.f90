module particules

use iso_c_binding
use zone
use quietstart

implicit none

integer, private :: ipart 
integer, private :: i, j

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine interpol_eb( tm1, ele )

type(particle) :: ele
type(mesh_fields) :: tm1
real(c_double) :: a1, a2, a3, a4
real(c_double) :: xp, yp, dxp, dyp

do ipart=1,nbpart

   xp = ele%pos(ipart,1) / dx
   yp = ele%pos(ipart,2) / dy

   i = floor( xp )
   j = floor( yp )

   dxp = xp - i
   dyp = yp - j

   a1 = (1-dxp) * (1-dyp) 
   a2 = dxp * (1-dyp) 
   a3 = dxp * dyp 
   a4 = (1-dxp) * dyp 

   ele%epx(ipart) = a1 * tm1%ex(i,j) + a2 * tm1%ex(i+1,j) &
        & + a3 * tm1%ex(i+1,j+1) + a4 * tm1%ex(i,j+1) 
   ele%epy(ipart) = a1 * tm1%ey(i,j) + a2 * tm1%ey(i+1,j) &
        & + a3 * tm1%ey(i+1,j+1) + a4 * tm1%ey(i,j+1) 
   ele%bpz(ipart) =  a1 * tm1%bz(i,j) + a2 * tm1%bz(i+1,j) &
        & + a3 * tm1%bz(i+1,j+1) + a4 * tm1%bz(i,j+1) 
end do

end subroutine interpol_eb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_vitesse( ele )

type (particle) :: ele
real(c_double) :: dum
real(c_double) :: tantheta, sintheta
real(8) :: hdt, v1, v2, e1, e2, b3

hdt = 0.5 * dt

do ipart = 1, nbpart

   v1 = ele%vit(ipart,1)
   v2 = ele%vit(ipart,2)
   e1 = ele%epx(ipart)
   e2 = ele%epx(ipart)
   b3 = ele%bpz(ipart)

   v1 = v1 + hdt * e1
   v2 = v2 + hdt * e2

   tantheta = hdt * b3
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   v1 = v1 + v2 * tantheta
   v2 = v2 - v1 * sintheta
   v1 = v1 + v2 * tantheta

   ele%vit(ipart,1) = v1 + hdt * e1
   ele%vit(ipart,2) = v2 + hdt * e2

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
real(c_double) :: coef

ele%pos = ele%pos + ele%vit * dt * coef

do ipart=1,nbpart
   ele%pos(ipart,1) = modulo(ele%pos(ipart,1), dimx)
   ele%pos(ipart,2) = modulo(ele%pos(ipart,2), dimy)
end do   

end subroutine avancee_part

subroutine calcul_j_cic( ele, tm, tm1 )

type(particle) :: ele
type(mesh_fields) :: tm, tm1
real(c_double) :: a1, a2, a3, a4, dum, xp, yp
integer :: ip1, jp1

tm1%jx = 0.d0
tm1%jy = 0.d0

do ipart=1,nbpart
   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)
   i = floor( xp / dx )
   j = floor( yp / dy )
   dum = dimx * dimy / (dx*dy) / nbpart
   a1 = (x(i+1)-xp) * (y(j+1)-yp) * dum
   a2 = (xp-x(i)) * (y(j+1)-yp) * dum
   a3 = (xp-x(i)) * (yp-y(j)) * dum
   a4 = (x(i+1)-xp) * (yp-y(j)) * dum
   dum = ele%vit(ipart,1) / (dx*dy) !charge unite = 1

   ip1 = mod(i+1,nx)
   jp1 = mod(j+1,ny)

   tm1%jx(i,j)     = tm1%jx(i,j)     + a1*dum  
   tm1%jx(ip1,j)   = tm1%jx(ip1,j)   + a2*dum 
   tm1%jx(ip1,jp1) = tm1%jx(ip1,jp1) + a3*dum 
   tm1%jx(i,jp1)   = tm1%jx(i,jp1)   + a4*dum 
   dum = ele%vit(ipart,2) / (dx*dy) 
   tm1%jy(i,j)     = tm1%jy(i,j)     + a1*dum  
   tm1%jy(ip1,j)   = tm1%jy(ip1,j)   + a2*dum 
   tm1%jy(ip1,j+1) = tm1%jy(ip1,jp1) + a3*dum 
   tm1%jy(i,jp1)   = tm1%jy(i,jp1)   + a4*dum 
end do

do i=0,nx
   tm1%jx(i,ny) = tm1%jx(i,0)
   tm1%jy(i,ny) = tm1%jy(i,0)
end do
do j=0,ny
   tm1%jx(nx,j) = tm1%jx(0,j)
   tm1%jy(nx,j) = tm1%jy(0,j)
end do


do i=0,nx-1
do j=0,ny
   tm%jx(i,j) = 0.5 * (tm1%jx(i,j)+tm1%jx(i+1,j))
end do
end do

do i=0,nx
do j=0,ny-1
   tm%jy(i,j) = 0.5 * (tm1%jy(i,j)+tm1%jy(i,j+1))
end do
end do

end subroutine calcul_j_cic

subroutine plasma( ele )

type (particle) :: ele
real(c_double) :: speed, theta, vth, n
real(c_double) :: a, b, eps, R
integer :: k

eps = 1.d-12

vth =  1.
nbpart = 100*(nx)*(ny)
n = 1.d0/nbpart

allocate(ele%pos(nbpart,2))
allocate(ele%vit(nbpart,2))
allocate(ele%epx(nbpart))
allocate(ele%epy(nbpart))
allocate(ele%bpz(nbpart))

do k=0,nbpart-1

   speed = vth * sqrt(-2 * log( (k+0.5)*n ))

   theta = trinary_reversing( k ) * 2 * pi

   a = 0; b = dimx ! 2*pi/kx 
   R = bit_reversing( k )
   call dichotomie_x(a,b,R,eps) 
   ele%pos(k+1,1) = a
   ele%pos(k+1,2) = dimy * penta_reversing( k ) 

   ele%vit(k+1,1) = speed * cos(theta)  !
   ele%vit(k+1,2) = speed * sin(theta)  !

enddo

end subroutine plasma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module particules
