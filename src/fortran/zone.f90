module zone

integer, parameter :: prec=8

type mesh_fields
   real(8), dimension(:,:), pointer :: ex, ey, bz, jx, jy
end type mesh_fields

type particle
   real(8), pointer :: pos(:,:)
   integer, pointer :: case(:,:)
   real(8), pointer :: vit(:,:)
   real(8), pointer :: epx(:)
   real(8), pointer :: epy(:)
   real(8), pointer :: bpz(:)
end type particle

real(8) :: pi 

integer :: nx, ny
integer :: nstep
integer :: nbpart

integer, private :: i, j

real(8) :: dt, alpha, kx
real(8) :: dx, dy
real(8), allocatable :: x(:)
real(8), allocatable :: y(:)
real(8) :: dimx, dimy
real(8) :: tfinal


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine init( )

implicit none

pi = 4. * atan(1.)
alpha = 0.1
kx = 0.5
dimx = 2*pi/kx
dimy = 1.0
nx = 128
ny = 16

allocate(x(0:nx+1))
allocate(y(0:ny+1))

dx = dimx / nx
dy = dimy / ny

do i=0,nx+1
   x(i) = i*dx 
enddo
do j=0,ny+1
   y(j) = j*dy
enddo

dt    = 0.01
nstep = 250

write(*,*) " dx = ", dx, " dy = ", dy, " dt = ", dt
write(*,*) " Nombre d'iteration nstep = ", nstep

end subroutine init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module zone
