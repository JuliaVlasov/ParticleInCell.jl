module zone

use iso_c_binding

implicit none

real(c_double) :: pi 

integer :: nx, ny
integer :: nbpart

integer, private :: i, j

real(c_double), allocatable :: x(:)
real(c_double), allocatable :: y(:)

real(c_double) :: dt, alpha, kx
real(c_double) :: dx, dy
real(c_double) :: dimx, dimy


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function mod1( x, y)
integer :: x
integer :: y

mod1 = modulo(x-1, y) + 1

end function mod1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init( )

implicit none

pi = 4. * atan(1d0)
alpha = 0.1
kx = 0.5
dimx = 2*pi/kx
dimy = 1.0
nx = 128
ny = 16

allocate(x(1:nx+1))
allocate(y(1:ny+1))

dx = dimx / nx
dy = dimy / ny

do i=1,nx+1
   x(i) = (i-1)*dx 
enddo
do j=1,ny+1
   y(j) = (j-1)*dy
enddo

dt    = 0.01

write(*,*) " dx = ", dx, " dy = ", dy, " dt = ", dt

end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
