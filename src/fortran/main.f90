program landau_damping

use pic

implicit none

integer, parameter :: nstep = 250
real(c_double) :: time(nstep+1)
real(c_double) :: energy(nstep+1)

integer :: istep

call vm2d2v( nstep, time, energy )

open(10,file='modeE.dat')
do istep = 1, nstep+1
    write(10,*) time(istep), energy(istep)
end do
close(10)
        

end program landau_damping
