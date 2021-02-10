$L_x,L_y$ domain dimensions and M,N are integers.
$
\omega = \sqrt{(\frac{M\pi}{L_x})^2+(\frac{N\pi}{L_y})^2}
$
$
B_z(x,y,t) =   - \cos(M \pi \frac{x}{L_x})  \cos(N \pi \frac{y    }{L_y}) \cos(\omega t)
$
$
E_x(x,y,t) = \frac{c^2 N \pi }{\omega Ly} cos(M \pi \frac{x}{L    _x}) \sin(N \pi  \frac{y}{L_y}) \sin(\omega t)
$
$
E_y(x,y,t) = - \frac{c^2 M \pi }{\omega Lx} \sin (M \pi \frac{    x}{L_x}) \cos (N \pi  \frac{y}{L_y}) \sin(\omega t)
$



```julia
using Pkg
Pkg.activate("/Users/navaro/JuliaProjects/ParticleInCell.jl")
```

```julia
using Plots
```

```julia
using Revise
```

```julia
using ParticleInCell
using LinearAlgebra
```

```julia
dimx, dimy = 2π, 2π
nx, ny = 128, 128
dt = 1e-3
nstep = 10
mesh = Mesh(dimx, nx, dimy, ny)
fdtd = FDTD(mesh)
ω = sqrt(2)

ex = zeros(nx, ny)
ey = zeros(nx, ny)
bz = zeros(nx, ny)
jx = zeros(nx, ny)
jy = zeros(nx, ny)

# Ex and Ey are set at t = 0.0
# Bz is set at  t = -dt/2

x = mesh.x[1:nx]
y = transpose(mesh.y[1:ny])
xc = ( x[1:nx-1] .+ x[2:nx] ) ./ 2
yc = transpose( y[1:ny-1] .+ y[2:ny] ) ./ 2

t = 0
ex .= 0
ey .= 0
t = 0.5dt
bz .= -cos.(x) .* cos.(y) .* cos(ω * t)
fdtd.bz .= cos.(xc) .* cos.(yc) .* cos(ω * t)

sol_ex = cos.(x) .* sin.(y) ./ ω
sol_ey = -sin.(x) .* cos.(y) ./ ω
sol_bz = -cos.(x) .* cos.(y)
 
ampere_maxwell!(ex, ey, fdtd, mesh, dt)
t = dt
@show maximum(abs.(ex .- sol_ex .* sin(ω * t))) 
@show maximum(abs.(ey .- sol_ey .* sin(ω * t))) 

faraday!(bz, maxwell, mesh, dt)
t = 1.5dt

@show maximum(abs.(bz .- sol_bz .* cos(ω * t)))
surface(bz)
```

```julia

```
