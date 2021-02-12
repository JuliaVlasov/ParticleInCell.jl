program VM_2D

use zone
use particules
use maxwell

implicit none

type(mesh_fields) :: f0, f1
type(particle) :: p

real(kind=prec) :: time
integer :: istep
integer :: i, j

call init( )

! staggered grid for FDTD maxwell scheme

allocate(f0%ex(0:nx-1,0:ny))
allocate(f0%ey(0:nx,0:ny-1))
allocate(f0%bz(0:nx-1,0:ny-1))
allocate(f0%jx(0:nx-1,0:ny))
allocate(f0%jy(0:nx,0:ny-1))


! Regular grid for particle interpolation and deposition

allocate(f1%ex(0:nx,0:ny)) 
allocate(f1%ey(0:nx,0:ny))
allocate(f1%bz(0:nx,0:ny))
allocate(f1%jx(0:nx,0:ny))
allocate(f1%jy(0:nx,0:ny))

time  = 0.d0

istep = 1

f0%ex = 0.d0; f0%ey = 0.d0; f0%bz = 0.d0
do i=0,nx-1
   do j=0,ny
      f0%ex(i,j) = alpha/kx * sin(kx*x(i))
   enddo
enddo

call plasma( p ) 

do istep = 1, nstep

   if (istep > 1) call faraday( f0%ex, f0%ey, f0%bz, 0.5*dt )

   call decalage( f0, f1 )
   call interpol_eb( f1, p )

   call avancee_vitesse( p )

   call avancee_part( p, 0.5d0 )  ! x(n) --> x(n+1/2)
   call calcul_j_cic( p, f0, f1 )
   call avancee_part( p, 0.5d0 )  ! x(n+1/2) -- x(n+1)
        
   call faraday( f0%ex, f0%ey, f0%bz, 0.5*dt )
   call ampere( f0%ex, f0%ey, f0%bz, f0%jx, f0%jy, dt ) 

   time = time + dt
   print*,'time = ',time, ' nbpart = ', nbpart

   call modeE( f0, istep, time )

end do


contains

subroutine modeE( f, iplot, time )

type(mesh_fields) :: f
real(kind=prec) :: time, aux
integer :: iplot, i, j


open(34,file='modeE.dat',position="append")
if (iplot==1) rewind(34)
write(34,*) time, 0.5*log(sum(f%ex*f%ex)*dx*dy)
close(34)

end subroutine modeE

end program VM_2D
