module pic

use iso_c_binding, only: c_int64_t, c_double

implicit none

public: interpolation, deposition

integer, private :: ipart 
integer, private :: i, j


contains

subroutine interpolation( pdata, fdata, nx, dx, ny, dy ) bind(C, name="interpolation")

real(c_double) :: a1, a2, a3, a4, s
real(c_double) :: xp, yp

s = dx * dy

do ipart=1,nbpart

   xp = pdata(1, ipart) / dx
   yp = pdata(2, ipart) / dy

   i = floor( xp ) + 1
   j = floor( yp ) + 1

   dxp = i - 1 - xp
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

subroutine deposition( ele, tm, tm1 ) bind(C, name="deposition")

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

end subroutine deposition

end module pic
