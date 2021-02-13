module pic

    use iso_c_binding
    use zone
    use particules
    use maxwell

    implicit none

    public :: vm2d2v

contains

    subroutine vm2d2v(nstep, time, energy) bind( C, name="vm2d2v" )

        integer(c_int32_t), intent(in) :: nstep
        real(c_double), intent(out) :: time(:)
        real(c_double), intent(out) :: energy(:)

        real(c_double), allocatable :: f0(:,:,:), f1(:,:,:)
        real(c_double), allocatable :: j0(:,:,:), j1(:,:,:)
        real(c_double), allocatable :: p(:,:)
        
        integer :: istep
        integer :: i, j
        
        write(*,*) " Nombre d'iteration nstep = ", nstep
        call init( )
        
        ! staggered grid for FDTD maxwell scheme
        allocate(f0(3,1:nx,1:ny))
        allocate(j0(2,1:nx,1:ny))
        
        ! Regular grid for particle interpolation and deposition
        allocate(f1(3,1:nx,1:ny)) 
        allocate(j1(2,1:nx,1:ny)) 
        
        istep = 1
        
        f0 = 0d0
        f1 = 0d0
        do i=1,nx
           do j=1,ny+1
              f0(1,i,j) = alpha/kx * sin(kx*x(i))
           end do
        end do
        
        call plasma( p ) 
        
        time(1)  = 0.d0
        energy(1)  = compute_energy( f0 )
        
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
        
           time(istep+1) = time(istep) + dt
           energy(istep+1) =  compute_energy( f0 )
           print*, istep, time(istep+1), energy(istep+1)
        
        end do
        
    end subroutine vm2d2v

end module pic
