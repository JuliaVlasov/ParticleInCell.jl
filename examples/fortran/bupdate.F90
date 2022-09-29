program uapic_2d

    use mesh_fields_m
    use particles_m
    use poisson_2d_m
    use m6_interpolation_m
    use m6_compute_rho_m
    use ua_steps_m

    implicit none


    real(8),    parameter :: alpha    = 0.05d0 
    real(8),    parameter :: kx       = 0.5d0
    real(8),    parameter :: ky       = 1d0
    integer,    parameter :: nx       = 128
    integer,    parameter :: ny       = 64 
    real(8),    parameter :: eps      = 0.1d0
    integer(8), parameter :: nbpart   = 204800_8
    integer,    parameter :: ntau     = 16
   
    real(8)    :: t = 0d0
    complex(8) :: elt

    real(8) :: tfinal   = 1.0d0 

    type(mesh_t)      :: mesh
    type(fields_2d_t) :: fields
    type(particles_t) :: particles
    type(poisson_t)   :: poisson
    type(ua_t)        :: ua

    complex(8), allocatable :: xt(:,:,:)
    complex(8), allocatable :: xf(:,:,:)
    complex(8), allocatable :: yt(:,:,:)
    complex(8), allocatable :: yf(:,:,:)
    complex(8), allocatable :: fx(:,:,:)
    complex(8), allocatable :: fy(:,:,:)
    complex(8), allocatable :: gx(:,:,:)
    complex(8), allocatable :: gy(:,:,:)

    real(8), allocatable :: et(:,:,:)

    real(8) :: dimx
    real(8) :: dimy
    integer :: istep
    integer :: nstep
    integer :: m
    integer :: n

    complex(8) :: px
    complex(8) :: py
    real(8) :: dx
    real(8) :: dy
    real(8) :: dt

    pi    = 4d0 * atan(1d0)
    dimx  = 2d0*pi/kx
    dimy  = 2d0*pi/ky 

    call init_mesh( mesh, 0d0, dimx, nx, 0d0, dimy, ny )

    dx = mesh%dx
    dy = mesh%dy

    dt = pi / 2d0 / (2d0**3) 
    tfinal = pi / 2d0

    nstep  = floor(tfinal/dt)

    call init_fields(fields, mesh )
    
    call init_particles_2d( particles, nbpart, mesh, alpha, kx )

    call init_poisson( poisson, mesh )

    call init_ua( ua, ntau, eps, nbpart )

    allocate(et(ntau, 2, nbpart))
    allocate(xt(ntau, 2, nbpart))
    allocate(xf(ntau, 2, nbpart))
    allocate(yt(ntau, 2, nbpart))
    allocate(yf(ntau, 2, nbpart))
    allocate(fx(ntau, 2, nbpart))
    allocate(fy(ntau, 2, nbpart))
    allocate(gx(ntau, 2, nbpart))
    allocate(gy(ntau, 2, nbpart))

    call compute_rho_m6_real( fields, particles )

    call solve_poisson( poisson, fields )

    call interpolate_eb_m6_real( particles, fields )

    do istep = 1,nstep

        call preparation( ua, dt, particles, xt, yt) 

        call interpolation( particles, et, fields, ua, xt)

        call compute_f( fx, fy, ua, particles, xt, yt, et )

        call ua_step1( xt, xf, ua, particles, fx )
        call ua_step1( yt, yf, ua, particles, fy )

        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields ) 

        call interpolation( particles, et, fields, ua, xt)
        
        call compute_f( gx, gy, ua, particles, xt, yt, et )

        call ua_step2( xt, xf, ua, particles, fx, gx ) 
        call ua_step2( yt, yf, ua, particles, fy, gy )

        call deposition( particles, fields, ua, xt)

        call solve_poisson( poisson, fields )

        call interpolation( particles, et, fields, ua, xt)

        call compute_v( ua, particles, yt, yf )

        print*, " vx, vy : ", sum(particles%v(1,:)), sum(particles%v(2,:))


    end do


end program uapic_2d


