module pic

    use iso_c_binding
    implicit none

contains

pure integer function mod1( x, y)

integer(c_int32_t), intent(in) :: x
integer(c_int32_t), intent(in) :: y

mod1 = modulo(x-1, y) + 1

end function mod1

pure subroutine interpolation( nbpart, nx, ny, dx, dy, ex, ey, bz, p ) bind(C, name="interpolation")

integer(c_int32_t), intent(in) :: nbpart, nx, ny
real(c_double), intent(in) :: dx, dy
real(c_double), intent(inout) :: p(7,nbpart)
real(c_double), intent(in) :: ex(nx+1,ny+1)
real(c_double), intent(in) :: ey(nx+1,ny+1)
real(c_double), intent(in) :: bz(nx+1,ny+1)

real(c_double) :: a1, a2, a3, a4
real(c_double) :: xp, yp, dxp, dyp
integer(c_int32_t) :: i, j, ipart

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

   p(5,ipart) = a1 * ex(i,j) + a2 * ex(i+1,j) + a3 * ex(i+1,j+1) + a4 * ex(i,j+1)
   p(6,ipart) = a1 * ey(i,j) + a2 * ey(i+1,j) + a3 * ey(i+1,j+1) + a4 * ey(i,j+1)
   p(7,ipart) = a1 * bz(i,j) + a2 * bz(i+1,j) + a3 * bz(i+1,j+1) + a4 * bz(i,j+1)

end do

end subroutine interpolation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine push_v( nbpart, p, dt ) bind(C, name="push_v")

integer(c_int32_t), intent(in) :: nbpart
real(c_double), intent(inout) :: p(7,nbpart)
real(c_double), intent(in) :: dt

real(c_double) :: tantheta, sintheta
real(c_double) :: hdt, v1, v2, e1, e2, b3
integer :: ipart

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

pure subroutine push_x( nbpart, dimx, dimy, p, dt ) bind(C, name="push_x")

integer(c_int32_t), intent(in) :: nbpart
real(c_double), intent(in) :: dimx, dimy, dt
real(c_double), intent(inout) :: p(7,nbpart)
real(c_double) ::  p1, p2
integer :: ipart

do ipart=1,nbpart
   p1 = p(1,ipart) + p(3,ipart) * dt 
   p2 = p(2,ipart) + p(4,ipart) * dt 
   if ( p1 > dimx ) p1 = p1 - dimx
   if ( p1 < 0.0  ) p1 = p1 + dimx
   if ( p2 > dimy ) p2 = p2 - dimy
   if ( p2 < 0.0  ) p2 = p2 + dimy
   p(1,ipart) = p1
   p(2,ipart) = p2
end do   

end subroutine push_x

pure subroutine deposition( nbpart, nx, ny, dx, dy, p, jx, jy ) bind(C, name="deposition")

integer(c_int32_t), intent(in) :: nbpart, nx, ny
real(c_double), intent(in)  :: dx, dy
real(c_double), intent(in)  :: p(7,nbpart)
real(c_double), intent(out) :: jx(nx+1,ny+1)
real(c_double), intent(out) :: jy(nx+1,ny+1)

real(c_double) :: a1, a2, a3, a4, w1, w2, xp, yp, dxp, dyp, factor
integer(c_int32_t) :: ipart, i, j

jx = 0.0
jy = 0.0

factor = real(nx * ny) / nbpart

do ipart=1,nbpart

   xp = p(1,ipart) / dx
   yp = p(2,ipart) / dy

   i = floor(xp) + 1
   j = floor(yp) + 1

   dxp = xp - i + 1
   dyp = yp - j + 1

   a1 = (1-dxp) * (1-dyp) 
   a2 = dxp * (1-dyp) 
   a3 = dxp * dyp 
   a4 = (1-dxp) * dyp 

   w1 = p(3,ipart) * factor
   w2 = p(4,ipart) * factor

   jx(i,j)     = jx(i,j)     + a1*w1  
   jx(i+1,j)   = jx(i+1,j)   + a2*w1 
   jx(i+1,j+1) = jx(i+1,j+1) + a3*w1 
   jx(i,j+1)   = jx(i,j+1)   + a4*w1 

   jy(i,j)     = jy(i,j)     + a1*w2  
   jy(i+1,j)   = jy(i+1,j)   + a2*w2 
   jy(i+1,j+1) = jy(i+1,j+1) + a3*w2 
   jy(i,j+1)   = jy(i,j+1)   + a4*w2 

end do

end subroutine deposition

end module pic
