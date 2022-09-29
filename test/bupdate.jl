using Test
using FFTW
using LinearAlgebra

include("../src/meshfields.jl")
include("../src/ua_type.jl")
include("../src/particles.jl")
include("../src/compute_rho.jl")
include("../src/gnuplot.jl")
include("../src/integrate.jl")
include("../src/interpolation.jl")
include("../src/plasma.jl")
include("../src/poisson.jl")
include("../src/read_particles.jl")
include("../src/ua_steps.jl")

function test_pic2d( ntau )

    kx       = 0.50
    ky       = 1.0
    dimx     = 2π/kx
    dimy     = 2π/ky 
    nx       = 128	
    ny       = 64 
    tfinal   = 1.0 

    t = 0.

    xmin, xmax = 0.0, dimx
    ymin, ymax = 0.0, dimy

    mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )

    dx, dy = mesh.dx, mesh.dy

    dt = π / 2 / (2^3) #Used to determine N in MRC, not the macro-step
    tfinal = π / 2

    nstep  = trunc(Int64, tfinal/dt)

    fields = MeshFields( mesh )
    
    #particles = read_particles( "particles.dat", mesh )
    particles = plasma( mesh, 204800 )

    nbpart = particles.nbpart

    poisson! = Poisson(mesh)

    ε = 0.1
    
    ua = UA( ntau, ε, nbpart )

    tau  = ua.tau
    ltau = ua.ltau

    et  = zeros(Float64, (ntau, 2, nbpart))

    xt  = zeros(ComplexF64, (ntau, 2, nbpart))
    x̃t  = zeros(ComplexF64, (ntau, 2, nbpart))
    yt  = zeros(ComplexF64, (ntau, 2, nbpart))
    ỹt  = zeros(ComplexF64, (ntau, 2, nbpart))

    ftau = plan_fft(xt,  1)

    fx = zeros(ComplexF64, (ntau, 2, nbpart))
    fy = zeros(ComplexF64, (ntau, 2, nbpart))

    gx = zeros(ComplexF64, (ntau, 2, nbpart))
    gy = zeros(ComplexF64, (ntau, 2, nbpart))

    compute_rho_m6!( fields, particles )

    nrj =  poisson!( fields )

    interpol_eb_m6!( particles, fields )

    for istep = 1:nstep

        preparation!( ua, dt, particles, xt, yt) 

        update_particles_e!( particles, et, fields, ua, xt)

        #  prediction

        compute_f!( fx, fy, ua, particles, xt, yt, et )

        mul!(x̃t, ftau, xt)
        ua_step!( xt, x̃t, ua, particles, fx )

        mul!(ỹt, ftau, yt)
        ua_step!( yt, ỹt, ua, particles, fy )

        ifft!(xt,1)
        ifft!(yt,1)

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields ) 

        update_particles_e!( particles, et, fields, ua, xt)

        # correction

        compute_f!( gx, gy, ua, particles, xt, yt, et )

        ua_step!( xt, x̃t, ua, particles, fx, gx ) 

        ua_step!( yt, ỹt, ua, particles, fy, gy )

        ifft!(xt,1)

        update_particles_x!( particles, fields, ua, xt)

        nrj = poisson!( fields )

        update_particles_e!( particles, et, fields, ua, xt)

        compute_v!( yt, particles, ua )

        @show @views sum(particles.v[1,:]), sum(particles.v[2,:])

    end


    true

end 

const ntau = 16

@time test_pic2d( ntau )
