program landau_damping

use pic

implicit none

integer, parameter :: nstep = 250
integer :: istep
real(c_double), allocatable :: time(:)
real(c_double), allocatable :: energy(:)

allocate(time(nstep+1))
allocate(energy(nstep+1))

call vm2d2v( nstep, time, energy )

open(10,file='modeE.dat',position="append")
do istep = 1, nstep+1
    write(10,*) time(istep), energy(istep)
end do
close(10)

end program landau_damping
