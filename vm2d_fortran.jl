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
    landau_sampling!( particles, nbpart, alpha, kx )

    fdtd = FDTD(mesh)
    for i=1:nx, j=1:ny
        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
        fdtd.ey[i,j] = 0.0
        fdtd.bz[i,j] = 0.0
    end
    time = 0
    energy = Float64[compute_energy(fdtd)]
    t = Float64[time]
    
    reset_timer!()

    @showprogress 1 for istep in 1:nstep
    
       if istep > 1
           @timeit to "fdtd" faraday!( fdtd, 0.5dt ) 
       end
       @timeit to "interpolation" interpolation!(particles, fdtd)
       @timeit to "pushv" f90_push_v!( particles, nbpart, dt )
       @timeit to "pushx" f90_push_x!( particles, nbpart, dimx, dimy, 0.5dt) 
       @timeit to "deposition" deposition!( fdtd, particles)
       @timeit to "pushx" f90_push_x!( particles, nbpart, dimx, dimy, 0.5dt) 
       @timeit to "fdtd" faraday!(fdtd, 0.5dt)
       @timeit to "fdtd" ampere_maxwell!(fdtd, dt)
       time = time + dt
       push!(t, time)
       push!(energy, compute_energy(fdtd))
    
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
