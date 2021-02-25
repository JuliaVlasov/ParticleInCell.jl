# Maxwell solver

```math
$L_x,L_y$ domain dimensions and M,N are integers.

$$
\omega = \sqrt{(\frac{M\pi}{L_x})^2+(\frac{N\pi}{L_y})^2}
$$

$$
B_z(x,y,t) =   - \cos(M \pi \frac{x}{L_x})  \cos(N \pi \frac{y}{L_y}) \cos(\omega t)
$$

$$
E_x(x,y,t) = \frac{c^2 N \pi }{\omega Ly} cos(M \pi \frac{x}{L_x}) \sin(N \pi  \frac{y}{L_y}) \sin(\omega t)
$$

$$
E_y(x,y,t) = - \frac{c^2 M \pi }{\omega Lx} \sin (M \pi \frac{x}{L_x}) \cos (N \pi  \frac{y}{L_y}) \sin(\omega t)
$$
```




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
mesh = TwoDGrid( dimx, nx, dimy, ny )
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

    x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])
    y = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose

    maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))
    
    
    @gif for istep = 1:nstep # Loop over time
    
        faraday!(maxwell, mesh, dt)     
    
        ampere_maxwell!(maxwell, mesh, dt) 
    
        surface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1))

    end every (nstep ÷ 100)
    
    
end

run(mesh, maxwell, 2000)
```
