module poisson_2d_m

use, intrinsic :: iso_c_binding 

use mesh_fields_m

implicit none

private

include 'fftw3.f03'

type, public :: poisson_t

    type(mesh_t)            :: mesh
    complex(8), allocatable :: kx(:,:)
    complex(8), allocatable :: ky(:,:)
    complex(8), allocatable :: rho(:,:)
    complex(8), allocatable :: ex(:,:)
    complex(8), allocatable :: ey(:,:)
    integer(8)              :: fw
    integer(8)              :: bw

end type poisson_t

public :: init_poisson, solve_poisson

contains

    subroutine init_poisson( self, mesh )
        
        type(poisson_t)      :: self
        type(mesh_t)         :: mesh

        complex(8)              :: kx0
        complex(8)              :: ky0
        real(8)                 :: pi
        real(8),    allocatable :: tmp(:,:)
        integer                 :: nx
        integer                 :: ny
        complex(8)              :: kx1
        complex(8), allocatable :: k2(:,:)
        integer                 :: ik
        integer                 :: jk

        pi = 4d0 * atan(1d0)

        kx0 = cmplx(2d0 * pi / (mesh%xmax - mesh%xmin), 0d0, kind=8)
        ky0 = cmplx(2d0 * pi / (mesh%ymax - mesh%ymin), 0d0, kind=8)

        nx = mesh%nx
        ny = mesh%ny

        allocate(self%kx(nx/2+1,ny))
        allocate(self%ky(nx/2+1,ny))
        allocate(k2(nx/2+1,ny))

        do ik=1,mesh%nx/2+1
           kx1 = cmplx(ik-1,0d0,kind=8)*kx0
           do jk = 1,mesh%ny/2
              self%kx(ik,jk) = kx1
              self%ky(ik,jk) = cmplx(jk-1,0d0,kind=8)*ky0
           end do
           do jk = mesh%ny/2+1,mesh%ny
              self%kx(ik,jk) = kx1
              self%ky(ik,jk) = cmplx(jk-1-ny,0d0,kind=8)*ky0
           end do
        end do

        self%kx(1,1) = (1d0,0d0)
        k2 = self%kx * self%kx + self%ky * self%ky
        self%kx = self%kx / k2
        self%ky = self%ky / k2

        allocate(self%rho(nx/2+1,ny))
        allocate(self%ex(nx/2+1,ny))
        allocate(self%ey(nx/2+1,ny))
        allocate(tmp(mesh%nx, mesh%ny))

        call dfftw_plan_dft_r2c_2d(self%fw,nx,ny,tmp,self%rho,FFTW_ESTIMATE)
        call dfftw_plan_dft_c2r_2d(self%bw,nx,ny,self%rho,tmp,FFTW_ESTIMATE)

    end subroutine init_poisson

    subroutine solve_poisson ( self, fields )

        type(poisson_t)     :: self
        type(fields_2d_t)   :: fields

        integer               :: nx
        integer               :: ny

        nx = fields%mesh%nx
        ny = fields%mesh%ny

        call dfftw_execute_dft_r2c( self%fw, fields%rho(1:nx,1:ny), self%rho )

        self%ex = cmplx(0d0,-1d0) * self%kx * self%rho
        self%ey = cmplx(0d0,-1d0) * self%ky * self%rho

        call dfftw_execute_dft_c2r( self%bw, self%ex, fields%e(1,1:nx,1:ny) )
        call dfftw_execute_dft_c2r( self%bw, self%ey, fields%e(2,1:nx,1:ny) )

        fields%e(1,nx+1,:) = fields%e(1,1,:)
        fields%e(1,:,ny+1) = fields%e(1,:,1)
        fields%e(2,nx+1,:) = fields%e(2,1,:)
        fields%e(2,:,ny+1) = fields%e(2,:,1)

        fields%e = fields%e / real( nx * ny, kind = 8 )

    end subroutine solve_poisson

end module poisson_2d_m
