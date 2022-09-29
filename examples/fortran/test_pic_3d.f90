program test_pic_3d

use mesh_fields_m
use poisson_3d_m
use output_m
use particles_m
use compute_rho_cic_m
use interpolation_cic_m

implicit none

type(mesh_t)       :: mesh
type(fields_3d_t)  :: fields
type(poisson_3d_t) :: poisson
type(particles_t)  :: particles

integer,    parameter :: nx = 64, ny = 64, nz = 4
integer(8), parameter :: nbpart = nx * ny * 100       
real(8)               :: xmin, xmax, ymin, ymax, zmin, zmax
real(8)               :: pi
integer               :: i, j, k, m
real(8)               :: x, y, z
real(8)               :: dx, dy, dz
real(8)               :: err_x, err_y, err_z


pi = 4d0 * atan(1d0)

xmin = 0d0
xmax = 18d0
ymin = 0d0
ymax = 18d0
zmin = 0d0
zmax = 1d0

call init_mesh( mesh, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz )

dx = mesh%dx
dy = mesh%dy
dz = mesh%dz

call init_fields( fields, mesh )

call init_poisson( poisson, mesh)

print*, " number of particles = ", nbpart

call init_particles_3d( particles, nbpart, mesh )

call compute_rho_cic( fields, particles)

print*, sum(fields%rho) * dx * dy * dz, nbpart * particles%w

call solve_poisson( poisson, fields)

call write_data( 1, fields)

do k = 1, nz+1
    z = (k-1) * dz
    do j = 1, ny+1
        y = (j-1) * dy
        do i = 1, nx+1
            x = (i-1) * dx
            fields%e(1,i,j,k) = sin(x)
            fields%e(2,i,j,k) = sin(y)
            fields%e(3,i,j,k) = sin(z)
        end do
    end do
end do

call interpolate_eb_cic( particles, fields)

err_x = 0d0
err_y = 0d0
err_z = 0d0

do m = 1, nbpart
    x = (i-1) * mesh%dx
    err_x = err_x + abs(particles%e(1,m) - sin(particles%x(1,m)))
    err_y = err_y + abs(particles%e(2,m) - sin(particles%x(2,m)))
    err_z = err_z + abs(particles%e(3,m) - sin(particles%x(3,m)))
end do

print*, err_x/nbpart, err_y/nbpart, err_z/nbpart

end
