# Vlasov-Maxwell

## Landau damping

```@setup vm2d
using Plots, LinearAlgebra, ProgressMeter
```

```@example vm2d
using ParticleInCell

nx        = 128  # nombre de pts suivant x
ny        = 16   # nombre de pts suivant y

alpha = 0.1
kx = 0.5
ky = 0.
dimx = 2*pi/kx
dimy = 1  
poids = dimx * dimy 

mesh = Mesh( dimx, nx, dimy, ny)
fdtd = FDTD(mesh)

ex = zeros(nx,ny)
ey = zeros(nx,ny)
bz = zeros(nx,ny)
jx = zeros(nx,ny)
jy = zeros(nx,ny)

time  = 0

for i=1:nx, j=1:ny
    fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
end
surface!(mesh.x[1:nx], mesh.y[1:ny], fdtd.ex )
```

```@example vm2d
nbpart = 100*nx*ny
particles = Particles(nbpart)
landau_sampling!( particles, alpha, kx )
update_cells!( particles, mesh )

p = plot(layout=4)
histogram!(p[1], particles.pos[1,:], normalize=true, label="x")
histogram!(p[2], particles.pos[2,:], normalize=true, label="y")
histogram!(p[3], particles.vit[1,:], normalize=true, label="vx")
histogram!(p[4], particles.vit[2,:], normalize=true, label="vy")
```

```@example vm2d
compute_current!( jx, jy, fdtd, particles, mesh)

p = plot(layout=2)
surface!(p[1], jx )
surface!(p[2], jy)
```

```@example vm2d
function run( ex, ey, bz, jx, jy, particles, mesh, nstep, dt)
    
    alpha = 0.1
    kx = 0.5
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
    
    for istep in 1:nstep
    
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
```

```@example vm2d
dt = 0.01
nstep = 250
t, energy = run( ex, ey, bz, jx, jy, particles, mesh, nstep, dt)
plot(t, energy, m=:o)
```
