using LinearAlgebra, ProgressMeter
using TimerOutputs, Sobol

include("src/mesh.jl")
include("src/particles.jl")
include("src/fdtd.jl")
include("src/pic.jl")

const to = TimerOutput()

function run( nstep; npm :: Int = 100 )
    
    dt = 0.01
    alpha = 0.1
    kx   = 0.5
    dimx = 2*pi/kx
    dimy = 1.0  
    nx   = 128  # nombre de pts suivant x
    ny   = 16   # nombre de pts suivant y
    mesh = Mesh( dimx, nx, dimy, ny)
    dx, dy = mesh.dx, mesh.dy

    @show nbpart = npm*nx*ny
    particles = zeros(7,nbpart)
    landau_sampling!( particles, alpha, kx )

    fdtd = FDTD(mesh)
    for i=1:nx, j=1:ny+1
        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
    end
    time = 0
    energy = Float64[compute_energy(fdtd, mesh)]
    t = Float64[time]
    
    reset_timer!()

    @showprogress 1 for istep in 1:nstep
    
       if istep > 1
           @timeit to "fdtd" faraday!( fdtd, mesh, 0.5dt ) 
       end

       update_fields!(mesh, fdtd)
       @timeit to "interpolation" f90_interpolation!(particles, mesh)
       @timeit to "pushv" f90_push_v!( particles, dt )
       @timeit to "pushx" f90_push_x!( particles, mesh, 0.5dt) 
       @timeit to "deposition" f90_compute_current!( mesh, particles)
       @timeit to "pushx" f90_push_x!( particles, mesh, 0.5dt) 
       @timeit to "fdtd" faraday!(fdtd, mesh, 0.5dt)
       @timeit to "fdtd" ampere_maxwell!(fdtd, mesh, dt)
       time = time + dt
       push!(t, time)
       push!(energy, compute_energy(fdtd, mesh))
    
    end
   
    t, energy
    
end

t, energy = run( 1 )
@show nstep = 1000
@time t, energy = run( nstep, npm=500 )
show(to)

open("results_f90.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
