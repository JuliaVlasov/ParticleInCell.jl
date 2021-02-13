program vm2d2v

use zone
use particules
use maxwell

implicit none

real(c_double), allocatable :: f0(:,:,:), f1(:,:,:)
real(c_double), allocatable :: j0(:,:,:), j1(:,:,:)
real(c_double), allocatable :: p(:,:)
real(c_double) :: time

integer :: istep
integer :: i, j

call init( )

! staggered grid for FDTD maxwell scheme
allocate(f0(3,1:nx,1:ny))
allocate(j0(2,1:nx,1:ny))

! Regular grid for particle interpolation and deposition
allocate(f1(3,1:nx,1:ny)) 
allocate(j1(2,1:nx,1:ny)) 

time  = 0.d0

istep = 1

f0 = 0d0
f1 = 0d0
do i=1,nx
   do j=1,ny+1
      f0(1,i,j) = alpha/kx * sin(kx*x(i))
   end do
end do

call plasma( p ) 

call modeE( f0, 0, time )

do istep = 1, nstep

   if (istep > 1) call faraday( f0, 0.5*dt )

   call decalage( f0, f1 )
   call interpolation( f1, p )

   call push_v( p )

   call push_x( p, 0.5d0 )  ! x(n) --> x(n+1/2)
   call deposition( p, j0, j1 )
   call push_x( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
        
   call faraday( f0, 0.5*dt )
   call ampere( f0, j0, dt ) 

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   call modeE( f0, istep, time )

end do


contains

subroutine modeE( f, iplot, time )

real(c_double) :: f(:,:,:)
real(c_double) :: time, aux
integer :: iplot, i, j

open(34,file='modeE.dat',position="append")
if (iplot==0) rewind(34)
write(34,*) time, 0.5*log(sum(f(1,1:nx,1:ny)*f(1,1:nx,1:ny))*dx*dy)
close(34)

end subroutine modeE

end program vm2d2v
