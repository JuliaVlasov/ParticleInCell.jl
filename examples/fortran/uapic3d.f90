program ua_pic_3d

use mesh_fields_m
use poisson_3d_m
use output_m
use particles_m
use compute_rho_cic_m
use interpolation_cic_m

implicit none

type(mesh_t)       :: mesh
type(fields_3d_t)  :: f
type(poisson_3d_t) :: poisson
type(particles_t)  :: p

real(8), parameter :: delta = 3.d-3
real(8), parameter :: ep    = 0.5d0**10
real(8)            :: pi

integer               :: istep
integer               :: n, m
integer,    parameter :: nx = 64, ny = 64, nz = 4
integer(8), parameter :: nbpart = nx * ny * 100       
real(8)               :: xmin(3), xmax(3)
real(8)               :: dx, dy, dz
real(8)               :: err_x, err_y, err_z
integer               :: Nmrc
integer               :: N0mrc
integer               :: Nmrcm
integer               :: nstep
real(8)               :: Bm(3), Ee(3)
real(8)               :: dt
real(8)               :: time
real(8)               :: tfinal
real(8)               :: alpha
real(8)               :: beta

pi = 4d0 * atan(1d0)

xmin = [00d0, 00d0, 00d0]
xmax = [18d0, 18d0, 01d0]

tfinal = pi

Nmrc   = 2**7
Nmrcm  = 2**7
N0mrc  = nint(tfinal/ep/(2d0*pi)/Nmrc)
nstep  = 1

if (N0mrc == 0 .or. N0mrc == 1) then
    dt    = ep*(2d0*pi)/Nmrc
    nstep = nint(tfinal/dt)
else
    alpha = 0.5d0*(1d0+1d0/N0mrc)*ep*N0mrc
    beta  = 0.5d0*(1d0-1d0/N0mrc)*ep*N0mrc
    dt    = (2d0*pi)/Nmrcm
endif

call init_mesh( mesh, xmin(1), xmax(1), nx, &
                      xmin(2), xmax(2), ny, &
                      xmin(3), xmax(3), nz )

dx = mesh%dx
dy = mesh%dy
dz = mesh%dz

call init_fields( f, mesh )

call init_poisson( poisson, mesh)

print*, " number of particles = ", nbpart

call init_particles_3d( p, nbpart, mesh )

call compute_rho_cic( f, p)

print*, sum(f%rho) * dx * dy * dz, nbpart * p%w

call solve_poisson( poisson, f)

call write_data( 1, f)

call interpolate_eb_cic( p, f)

time  = 0.0_8

print"(a,9g10.3)", ' N0 = ', N0mrc, nstep

if (N0mrc == 0 .or. N0mrc == 1) then

    do istep = 1, nstep !*** Loop over time

        call push_particles(p, 0.5d0*dt)

        call compute_rho_cic(f, p)

        call solve_poisson(poisson, f)

        call interpolate_eb_cic( p, f )

        do m = 1, p%nbpart

            Bm(1) = (p%x(2,m)-9.d0)*delta/sqrt(1d0 &
                  +((p%x(1,m)-9d0)**2+(p%x(2,m)-9.d0)**2)*delta**2)

            Bm(2) =-(p%x(1,m)-9.d0)*delta/sqrt(1d0 &
                +((p%x(1,m)-9d0)**2+(p%x(2,m)-9.d0)**2)*delta**2)

            Bm(3) = 1d0/sqrt(1d0+((p%x(1,m)-9d0)**2+(p%x(2,m)-9d0)**2)*delta**2)

            Ee = p%e(:,m)

            p%v(:,m) = cos(dt/ep) * p%v(:,m) &
                     + sin(dt/ep) * cross(p%v(:,m),Bm) &
                     + ep * sin(dt/ep) * Ee &
                     + (dt-ep*sin(dt/ep)) * &
                       (Bm(1) * Ee(1)+Bm(2)*Ee(2) &
                     + Bm(3) * Ee(3)) * Bm &
                     + (ep - ep*cos(dt/ep))*cross(Ee,Bm) &
                     + (1.d0-cos(dt/ep))*(Bm(1)*p%v(1,m) &
                     + Bm(2)*p%v(2,m)+Bm(3)*p%v(3,m))*Bm
        end do

        p%x = p%x + 0.5d0 * dt * p%v

        time = time + dt

    end do ! next time step

else

    do istep = 1, Nmrc

        do n=1,Nmrcm

            call push_particles(p, 0.5d0*dt*alpha)

            call compute_rho_cic(f, p)
            call solve_poisson(poisson, f)
            call interpolate_eb_cic( p, f )

            do m=1,p%nbpart

                Bm(1) = (p%x(2,m)-9d0)*delta/sqrt(1d0+((p%x(1,m)-9d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)
                Bm(2) =-(p%x(1,m)-9d0)*delta/dsqrt(1d0+((p%x(1,m)-9.d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)
                Bm(3) = 1d0/dsqrt(1d0+((p%x(1,m)-9d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)

                Ee = p%e(:,m)

                p%v(:,m) = cos(dt)*p%v(:,m)+sin(dt)*cross(p%v(:,m),Bm) &
                          + alpha*sin(dt)*Ee &
                          + alpha*(dt-sin(dt)) &
                          * (Bm(1)*Ee(1)+Bm(2)*Ee(2)+Bm(3)*Ee(3))*Bm &
                          + alpha*(1d0-cos(dt))*cross(Ee,Bm) &
                          + (1d0-cos(dt))*(Bm(1)*p%v(1,m) &
                          +                Bm(2)*p%v(2,m) &
                          +                Bm(3)*p%v(3,m)) * Bm
            end do

            call push_particles(p, 0.5d0*dt*alpha)

        end do

        do n = 1, Nmrcm

            call push_particles(p, 0.5d0*dt*beta)

            call compute_rho_cic(f, p)
            call solve_poisson(poisson, f)
            call interpolate_eb_cic( p, f )

            do m=1,p%nbpart

                Bm(1) = (p%x(2,m)-9d0)*delta/sqrt(1d0+((p%x(m,1)-9d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)

                Bm(2) =-(p%x(1,m)-9d0)*delta/sqrt(1d0+((p%x(m,1)-9d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)

                Bm(3) = 1d0/sqrt(1d0+((p%x(1,m)-9d0)**2 &
                       +(p%x(2,m)-9d0)**2)*delta**2)

                Ee = p%e(:,m)

                p%v(:,m) =  cos( dt) * p%v(:,m) &
                          + sin(-dt) * cross(p%v(:,m),Bm) &
                          - beta * sin(-dt) * Ee &
                          - beta * (-dt-dsin(-dt)) &
                          * dot_product(Ee,Bm) * Bm &
                          - beta * (1d0-cos(dt)) * cross(Ee,Bm) &
                          + (1d0-cos(dt))*(Bm(1)*p%v(1,m) &
                          + Bm(2) * p%v(2,m) + Bm(3) * p%v(3,m))*Bm
            end do


            call push_particles(p, 0.5d0*dt*beta)

        end do

        time = time + N0mrc * (2d0*pi)

    end do

end if

contains

pure function cross(v,w) result(res)

    real(8), intent(in) :: v(3)
    real(8), intent(in) :: w(3)
    real(8)             :: res(3)

    res(1) = v(2) * w(3) - v(3) * w(2)
    res(2) = v(3) * w(1) - v(1) * w(3)
    res(3) = v(1) * w(2) - v(2) * w(1)

end function cross

subroutine push_particles( p, delta_t )

    type(particles_t)   :: p
    real(8), intent(in) :: delta_t

    integer :: m
    real(8) :: d(3)

    do m = 1, p%nbpart !periodic BC
    
        d = delta_t * p%v(:,m) - xmin(:)

        p%x(1,m) = xmin(1) + modulo(p%x(1,m) + d(1), xmax(1) - xmin(1))
        p%x(2,m) = xmin(2) + modulo(p%x(2,m) + d(2), xmax(2) - xmin(2))
        p%x(3,m) = xmin(3) + modulo(p%x(3,m) + d(3), xmax(3) - xmin(3))
        
    end do

end subroutine push_particles

end program ua_pic_3d
