using ProgressMeter
using TimerOutputs
using ParticleInCell

# +
const to = TimerOutput()

function run( nstep; npm = 100 )
    
    dt = 0.01
    alpha = 0.1
    kx   = 0.5
    dimx = 2pi/kx
    dimy = 1  
    nx   = 128  # nombre de pts suivant x
    ny   = 16   # nombre de pts suivant y
    mesh = TwoDGrid( dimx, nx, dimy, ny)
    dx, dy = mesh.dx, mesh.dy

    nbpart = npm*nx*ny
    println( " nbpart = $nbpart ")

    particles = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    sampler = LandauDamping(alpha, kx)

    sample!(particles, mesh, sampler)

    fdtd = FDTD(mesh)
    for i=1:nx, j=1:ny+1
        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
    end
    time = 0
    energy = Float64[compute_energy(fdtd, mesh)]
    t = Float64[time]

    kernel = CloudInCell()
    
    reset_timer!()
    @showprogress 1 for istep in 1:nstep
    
       if istep > 1
           @timeit to "fdtd" faraday!( fdtd, mesh, 0.5dt ) 
       end
       update_fields!(mesh, fdtd)
       #@timeit to "interpolation" interpolation!(particles, mesh)
       @timeit to "pushv" push_v!( particles, kernel, mesh, dt )
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "deposition" compute_current!( mesh, kernel, particles)
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "fdtd" faraday!(fdtd, mesh, 0.5dt)
       @timeit to "fdtd" ampere_maxwell!(fdtd, mesh, dt)
       time = time + dt
       push!(t, time)
       push!(energy, compute_energy(fdtd, mesh))
    
    end
   
    t, energy
    
end

# -

t, energy = run( 1 ) # trigger building
@show nstep = 1000
@time t, energy = run( nstep, npm = 500 )
show(to)

open("results_jl.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
