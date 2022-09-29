program test_poisson_3d

use mesh_fields_m
use poisson_3d_m
use output_m

implicit none

type(mesh_t)       :: mesh
type(fields_3d_t)  :: fields
type(poisson_3d_t) :: poisson

integer, parameter :: nx = 32, ny = 64, nz = 128
real(8)            :: xmin, xmax, ymin, ymax, zmin, zmax
real(8)            :: pi
integer            :: i, j, k
real(8)            :: x, y, z
real(8)            :: err_x, err_y, err_z

pi = 4d0 * atan(1d0)

xmin = 0d0
xmax = 2d0 * pi
ymin = 0d0
ymax = 4d0 * pi
zmin = 0d0
zmax = 6d0 * pi

call init_mesh( mesh, xmin, xmax, nx, ymin, ymax, ny, zmin, zmax, nz )

call init_fields( fields, mesh )

call init_poisson( poisson, mesh)

do k = 1, nz
    z = (k-1) * mesh%dz
    do j = 1, ny
        y = (j-1) * mesh%dy
        do i = 1, nx
            x = (i-1) * mesh%dx
            fields%rho(i,j,k) = - 3d0 * sin(x) * sin(y) * sin(z)
        end do
    end do
end do

call solve_poisson( poisson, fields)

call write_data( 1, fields)

err_x = 0d0
err_y = 0d0
err_z = 0d0
do k = 1, nz+1
    z = (k-1) * mesh%dz
    do j = 1, ny+1
        y = (j-1) * mesh%dy
        do i = 1, nx+1
            x = (i-1) * mesh%dx
            err_x = err_x + abs( fields%e(1,i,j,k) - cos(x)*sin(y)*sin(z))
            err_y = err_y + abs( fields%e(2,i,j,k) - cos(y)*sin(x)*sin(z))
            err_z = err_z + abs( fields%e(3,i,j,k) - cos(z)*sin(x)*sin(y))
        end do
    end do
end do

print*, err_x, err_y, err_z

end 
