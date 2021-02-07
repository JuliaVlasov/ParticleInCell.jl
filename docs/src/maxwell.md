# Maxwell solver


```@setup maxwell
using Plots
```

```@example maxwell
using ParticleInCell

dimx, dimy = 1, 1
nx, ny = 64, 64
e0 = 1
c = 1
md, nd = 2, 2  # number of modes along each direction
dt = 0.001
nstep = 1 ÷ dt
mesh = Mesh( dimx, nx, dimy, ny )
maxwell = MaxwellSolver( mesh, c, e0) 
omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)
ex = zeros(nx+1,ny+1)
ey = zeros(nx+1,ny+1)
bz = zeros(nx+1,ny+1)

x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])
y = transpose(0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]))

maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))
    
faraday!(maxwell, bz, dt)    
surface(bz, aspect_ratio=:equal, zlims=(-1,1))
```

- Ex and Ey are set at t = 0.0
- Bz is set at  t = -dt/2

```@example maxwell
function run(nstep)
    dimx, dimy = 1, 1
    nx, ny = 64, 64
    e0 = 1
    c = 1
    md, nd = 2, 2  # number of modes along each direction
    dt = 0.001
    mesh = Mesh( dimx, nx, dimy, ny )
    maxwell = MaxwellSolver( mesh, c, e0) 
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)
    ex = zeros(nx+1,ny+1)
    ey = zeros(nx+1,ny+1)
    bz = zeros(nx+1,ny+1)
    
    x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])
    y = transpose(0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]))

    maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))
    
    
    @gif for istep = 1:nstep # Loop over time
    
        faraday!(maxwell, bz, dt)     
    
        ampere_maxwell!(maxwell, ex, ey, dt) 
    
        #p = plot(layout=(1,2))
        #surface!( p[1,1], ex, aspect_ratio=:equal, zlims=(-1,1))
        #surface!( p[1,2], ey, aspect_ratio=:equal, zlims=(-1,1))
        surface(bz, aspect_ratio=:equal, zlims=(-1,1))

        end every (nstep ÷ 100)
    
    
end

run(2000)
```


