# Maxwell solver


```@setup maxwell
using Plots
```

```@example maxwell
using ParticleInCell

dimx, dimy = 1, 1
nx, ny = 64, 64
md, nd = 2, 2  
dt = 0.001
nstep = 1 ÷ dt
mesh = Mesh( dimx, nx, dimy, ny )
maxwell = FDTD( mesh ) 
omega = sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)

x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])
y = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose

maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))
    
surface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1))
```

- Ex and Ey are set at t = 0.0
- Bz is set at  t = -dt/2

```@example maxwell
function run(mesh, maxwell, nstep)

    ex = zeros(nx,ny)
    ey = zeros(nx,ny)
    bz = zeros(nx,ny)
    
    x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])
    y = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose

    maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))
    
    
    @gif for istep = 1:nstep # Loop over time
    
        faraday!(bz, maxwell, mesh, dt)     
    
        ampere_maxwell!(ex, ey, maxwell, mesh, dt) 
    
        surface(bz, aspect_ratio=:equal, zlims=(-1,1), clims=(-1,1))

    end every (nstep ÷ 100)
    
    
end

run(mesh, maxwell, 2000)
```


