using Plots, LinearAlgebra, ProgressMeter
using ParticleInCell

function run( nstep )
    
    dt = 0.01
    alpha = 0.1
    kx = 0.5
    dimx = 2*pi/kx
    dimy = 1  
    nx     = 128  # nombre de pts suivant x
    ny     = 16   # nombre de pts suivant y
    mesh = Mesh( dimx, nx, dimy, ny)
    nbpart = 100*nx*ny
    ex = zeros(nx,ny)
    ey = zeros(nx,ny)
    bz = zeros(nx,ny)
    jx = zeros(nx,ny)
    jy = zeros(nx,ny)
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
    energy = [0.5 * log( sum( fdtd.ex.^2) * mesh.dx * mesh.dy)]
    t = Float64[time]
    
    @showprogress 1 for istep in 1:nstep
    
       istep > 1 && faraday!( bz, fdtd, mesh, 0.5dt ) 
       interpol_eb!( ex, ey, bz, particles, mesh )
       push_v!( particles, dt )
       push_x!( particles, mesh, 0.5dt) 
       compute_current!( jx, jy, fdtd, particles, mesh)
       push_x!( particles, mesh, 0.5dt) 
       faraday!(bz, fdtd, mesh, 0.5dt)
       ampere_maxwell!(ex, ey, fdtd, mesh, dt)
       time = time + dt
       push!(t, time)
       push!(energy, 0.5 * log( sum(fdtd.ex.^2) * mesh.dx * mesh.dy))
    
    end
   
    t, energy
    
end

nstep = 1
t, energy = run( nstep )
nstep = 250
@time t, energy = run( nstep )
plot(t, energy, m=:o)
savefig("out.png")
