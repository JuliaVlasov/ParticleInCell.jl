module particules

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
real(kind = prec) :: a1, a2, a3, a4, s
real(kind = prec) :: xp, yp

s = dx * dy

do ipart=1,nbpart

   xp = ele%pos(ipart,1)
   yp = ele%pos(ipart,2)

   i = floor( xp / dx)
   j = floor( yp / dx)

   a1 = (x(i+1)-xp) * (y(j+1)-yp)  / s
   a2 = (xp-x(i)) * (y(j+1)-yp)  / s
   a3 = (xp-x(i)) * (yp-y(j)) / s
   a4 = (x(i+1)-xp) * (yp-y(j))/ s

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
real(kind=prec) :: dum
real(kind=prec) :: tantheta, sintheta

do ipart = 1, nbpart

   dum = 0.5 * dt
   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*ele%epx(ipart)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*ele%epy(ipart)

   tantheta = dum * ele%bpz(ipart)
   sintheta = 2.0 * tantheta / ( 1. + tantheta*tantheta)

   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta
   ele%vit(ipart,2) = ele%vit(ipart,2) - ele%vit(ipart,1)*sintheta
   ele%vit(ipart,1) = ele%vit(ipart,1) + ele%vit(ipart,2)*tantheta

   ele%vit(ipart,1) = ele%vit(ipart,1) + dum*ele%epx(ipart)
   ele%vit(ipart,2) = ele%vit(ipart,2) + dum*ele%epy(ipart)

end do

end subroutine avancee_vitesse

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine avancee_part( ele, coef )  !Avancee de coef * dt

type(particle) :: ele
real(kind=prec) :: coef

ele%pos = ele%pos + ele%vit * dt * coef

do ipart=1,nbpart
   ele%pos(ipart,1) = modulo(ele%pos(ipart,1), dimx)
   ele%pos(ipart,2) = modulo(ele%pos(ipart,2), dimy)
end do   

end subroutine avancee_part

subroutine calcul_j_cic( ele, tm, tm1 )

type(particle) :: ele
type(mesh_fields) :: tm, tm1
real(kind=prec) :: a1, a2, a3, a4, dum, xp, yp
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
real(kind=prec) :: speed, theta, vth, n
real(kind=prec) :: a, b, eps, R
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
