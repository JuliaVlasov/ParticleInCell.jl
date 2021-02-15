using LinearAlgebra, ProgressMeter
using ParticleInCell
using TimerOutputs

const to = TimerOutput()

function run( nstep )
    
    dt = 0.01
    alpha = 0.1
    kx   = 0.5
    dimx = 2*pi/kx
    dimy = 1  
    nx   = 128  # nombre de pts suivant x
    ny   = 16   # nombre de pts suivant y
    mesh = Mesh( dimx, nx, dimy, ny)
    dx, dy = mesh.dx, mesh.dy

    @show nbpart = 100*nx*ny
    particles = Particles(nbpart)
    landau_sampling!( particles, alpha, kx )

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
       @timeit to "interpolation" interpol_eb!(particles, fdtd)
       @timeit to "pushv" push_v!( particles, dt )
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "deposition" compute_current!( fdtd, particles)
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "fdtd" faraday!(fdtd, 0.5dt)
       @timeit to "fdtd" ampere_maxwell!(fdtd, dt)
       time = time + dt
       push!(t, time)
       push!(energy, compute_energy(fdtd))
    
    end
   
    t, energy
    
end

t, energy = run( 1 )
@show nstep = 250
@time t, energy = run( nstep )
show(to)

open("results_jl.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
