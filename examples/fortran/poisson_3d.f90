module poisson_3d_m

use mesh_fields_m

use, intrinsic :: iso_c_binding

implicit none

include 'fftw3.f03'


type, public :: poisson_3d_t

    type(mesh_t)            :: mesh
    complex(8), allocatable :: rhs(:,:,:)
    complex(8), allocatable :: psi(:,:,:)
    integer(8)              :: fw
    integer(8)              :: bw

end type poisson_3d_t


contains

subroutine init_poisson( poisson, mesh)

  type(poisson_3d_t)  :: poisson
  type(mesh_t)        :: mesh

  integer :: nx, ny, nz
  
  nx = mesh%nx
  ny = mesh%ny
  nz = mesh%nz

  allocate(poisson%rhs(nx,ny,nz))
  allocate(poisson%psi(nx,ny,nz))

  call dfftw_plan_dft_3d(poisson%fw, nx, ny, nz, poisson%rhs, poisson%psi, &
                         fftw_forward, fftw_estimate)

  call dfftw_plan_dft_3d(poisson%bw, nx, ny, nz, poisson%psi, poisson%rhs, &
                         fftw_backward, fftw_estimate)
  
end subroutine init_poisson

subroutine solve_poisson(poisson, fields)

  type(poisson_3d_t) :: poisson
  type(fields_3d_t)  :: fields

  integer  :: nx, ny, nz
  integer  :: ind_x, ind_y, ind_z
  real(8)  :: kx, ky, kz
  real(8)  :: dimx, dimy, dimz
  real(8)  :: pi
  integer  :: i, j, k

  pi = 4d0 * atan(1d0)

  nx = fields%mesh%nx
  ny = fields%mesh%ny
  nz = fields%mesh%nz

  dimx = fields%mesh%xmax - fields%mesh%xmin
  dimy = fields%mesh%ymax - fields%mesh%ymin
  dimz = fields%mesh%zmax - fields%mesh%zmin

  poisson%psi = cmplx(fields%rho(1:nx,1:ny,1:nz), 0d0, kind=8)

  call dfftw_execute_dft(poisson%fw, poisson%psi, poisson%rhs)
  
  do k=1,nz
      if (k<=nz/2) then
          ind_z = k-1
      else
          ind_z = -nz+(k-1)
      end if
      kz = 2d0 * pi * real(ind_z,8) / dimz
      do j=1,ny
          if (j<=ny/2) then
              ind_y = j-1
          else
              ind_y = -ny+(j-1)
          end if
          ky = 2d0 * pi * real(ind_y,8)  / dimy
          do i=1,nx
              if (i<=nx/2) then
                  ind_x = i-1
              else
                  ind_x = - nx+(i-1)
              end if
              kx = 2d0 * pi * real(ind_x,8)  / dimx
              if ( ind_x==0 .and. ind_y==0 .and. ind_z==0 ) then
                  poisson%rhs(i,j,k) = cmplx(0d0,0d0,kind=8)
              else
                  poisson%rhs(i,j,k) = -cmplx(0d0,kx,kind=8) * poisson%rhs(i,j,k) &
                          / cmplx(kx**2+ky**2+kz**2,0d0,kind=8)
              end if
         end do
      end do
  end do

  call dfftw_execute_dft( poisson%bw, poisson%rhs, poisson%psi )

  fields%e(1,1:nx,1:ny,1:nz) = real(poisson%psi)

  poisson%psi = cmplx(fields%rho(1:nx,1:ny,1:nz), 0d0, kind=8)

  call dfftw_execute_dft(poisson%fw, poisson%psi, poisson%rhs)

  do k=1,nz
      if (k<=nz/2) then
          ind_z = k-1
      else
          ind_z = -nz+(k-1)
      end if
      kz = 2d0 * pi * real(ind_z,8) / dimz
      do j=1,ny
          if (j<=ny/2) then
              ind_y = j-1
          else
              ind_y = -ny+(j-1)
          end if
          ky = 2d0 * pi * real(ind_y,8) / dimy
          do i=1,nx
              if (i<=nx/2) then
                  ind_x = i-1
              else
                  ind_x = - nx+(i-1)
              end if
              kx = 2d0 * pi * real(ind_x,8) / dimx
              if ( ind_x==0 .and. ind_y==0 .and. ind_z==0 ) then
                  poisson%rhs(i,j,k) = cmplx(0d0,0d0,kind=8)
              else
                  poisson%rhs(i,j,k) = -cmplx(0d0,ky,kind=8) * poisson%rhs(i,j,k) &
                                / cmplx(kx**2+ky**2+kz**2,0d0,kind=8)
              end if
          end do
    end do
  end do

  call dfftw_execute_dft( poisson%bw, poisson%rhs, poisson%psi )
  fields%e(2,1:nx,1:ny,1:nz) = real(poisson%psi, kind=8)

  poisson%psi = cmplx(fields%rho(1:nx,1:ny,1:nz), 0d0, kind=8)
  call dfftw_execute_dft(poisson%fw, poisson%psi, poisson%rhs)

  do k=1,nz
      if (k<=nz/2) then
          ind_z = k-1
      else
          ind_z = -nz+(k-1)
      endif
      kz = 2d0 * pi * real(ind_z,8) /dimz
      do j=1,ny
          if (j<=ny/2) then
              ind_y = j-1
          else
              ind_y = -ny+(j-1)
          endif
          ky = 2d0 * pi * real(ind_y,8) /dimy
          do i=1,nx
              if (i<=nx/2) then
                 ind_x = i-1
              else
                 ind_x = - nx+(i-1)
              end if
              kx = 2d0 * pi * real(ind_x,8) /dimx
              if ( ind_x==0 .and. ind_y==0 .and. ind_z==0 ) then
                  poisson%rhs(i,j,k) = cmplx(0d0,0d0,kind=8)
              else
                  poisson%rhs(i,j,k) = -cmplx(0d0,kz,kind=8) * poisson%rhs(i,j,k) &
                           / cmplx(kx**2+ky**2+kz**2,0d0,kind=8)
              end if
          end do
      end do
  enddo

  call dfftw_execute_dft( poisson%bw, poisson%rhs, poisson%psi )
  fields%e(3,1:nx,1:ny,1:nz) = real(poisson%psi, kind=8)

  !periodic boundary conditions

  fields%e(:,nx+1,:,:) = fields%e(:,1,:,:)
  fields%e(:,:,ny+1,:) = fields%e(:,:,1,:)
  fields%e(:,:,:,nz+1) = fields%e(:,:,:,1)

  fields%e = fields%e / real(nx*ny*nz,kind=8)

end subroutine solve_poisson

end module poisson_3d_m
