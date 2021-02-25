# Vlasov-Maxwell

## Landau damping

```@setup vm2d
using Plots, LinearAlgebra
```

```@example vm2d
using ParticleInCell

nx = 128 
ny = 16   

alpha = 0.1
kx = 0.5
ky = 0.
dimx = 2*pi/kx
dimy = 1  
poids = dimx * dimy 

mesh = TwoDGrid( dimx, nx, dimy, ny)
fdtd = FDTD(mesh)

time  = 0

for i=1:nx, j=1:ny
    fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
end
surface(fdtd.ex )
```

```@example vm2d
nbpart = 100*nx*ny
particles = zeros(7, nbpart)
landau_sampling!( particles, alpha, kx )

p = plot(layout=4)
histogram!(p[1], particles[1,:], normalize=true, label="x")
histogram!(p[2], particles[2,:], normalize=true, label="y")
histogram!(p[3], particles[3,:], normalize=true, label="vx")
histogram!(p[4], particles[4,:], normalize=true, label="vy")
```

```@example vm2d
compute_current!( mesh, particles)

p = plot(layout=2)
surface!(p[1], mesh.jx)
surface!(p[2], mesh.jy)
```

```@example vm2d
function run( fdtd, particles, mesh, nstep, dt)
    
    alpha = 0.1
    kx = 0.5
    landau_sampling!( particles, alpha, kx )
    for i=1:nx, j=1:ny
        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)
        fdtd.ey[i,j] = 0.0
        fdtd.bz[i,j] = 0.0
    end

    time = 0
    energy = Float64[0.5 * log( sum( fdtd.ex.^2) * mesh.dx * mesh.dy)]
    t = Float64[time]
    
    for istep in 1:nstep
    
       istep > 1 && faraday!( fdtd, mesh, 0.5dt ) 
       interpolation!( particles, mesh )
       push_v!( particles, dt )
       push_x!( particles, mesh, 0.5dt) 
       compute_current!( mesh, particles)
       push_x!( particles, mesh, 0.5dt) 
       faraday!(fdtd, mesh, 0.5dt)
       ampere_maxwell!(fdtd, mesh, dt)
       time = time + dt
       push!(t, time)
       push!(energy, 0.5 * log( sum(fdtd.ex.^2) * mesh.dx * mesh.dy))
    
    end
   
    t, energy
    
end
```

```@example vm2d
dt = 0.01
nstep = 1000
t, energy = run( fdtd, particles, mesh, nstep, dt)
plot(t, energy)
```
