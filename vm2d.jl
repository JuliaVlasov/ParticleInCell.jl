using Plots, LinearAlgebra, ProgressMeter
using ParticleInCell
using TimerOutputs

const to = TimerOutput()

function run( nstep )
    
    dt = 0.01
    alpha = 0.1
    kx = 0.5
    dimx = 2*pi/kx
    dimy = 1  
    nx     = 128  # nombre de pts suivant x
    ny     = 16   # nombre de pts suivant y
    mesh = Mesh( dimx, nx, dimy, ny)
    @show nbpart = 100*nx*ny
    jxy = zeros(2,nx,ny)
    particles = Particles(nbpart)
    landau_sampling!( particles, alpha, kx )
    update_cells!( particles, mesh )
    fdtd = FDTD(mesh)
    for i=1:nx, j=1:ny
        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
        fdtd.ey[i,j] = 0.0
        fdtd.bz[i,j] = 0.0
    end
    time = 0
    energy = Float64[0.5 * log( sum( fdtd.ex.^2) * mesh.dx * mesh.dy)]
    t = Float64[time]
    eb = zeros((3,nx,ny))
    
    reset_timer!()
    @showprogress 1 for istep in 1:nstep
    
       if istep > 1
           @timeit to "fdtd" faraday!( eb, fdtd, mesh, 0.5dt ) 
       end
       @timeit to "interpolation" interpol_eb!( eb, particles, mesh )
       @timeit to "pushv" push_v!( particles, dt )
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "deposition" compute_current!( jxy, fdtd, particles)
       @timeit to "pushx" push_x!( particles, mesh, 0.5dt) 
       @timeit to "fdtd" faraday!(eb, fdtd, 0.5dt)
       @timeit to "fdtd" ampere_maxwell!(eb, fdtd, dt)
       time = time + dt
       push!(t, time)
       push!(energy, 0.5 * log( sum(fdtd.ex.^2) * mesh.dx * mesh.dy))
    
    end

   
    t, energy
    
end

t, energy = run( 1 )
@show nstep = 250
@time t, energy = run( nstep )
plot(t, energy, m=:o)
savefig("out.png")
show(to)
