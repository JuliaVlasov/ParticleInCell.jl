# Landau damping


```@setup landau
using ParticleInCell
using Plots
using Random
```

## Test particle initialization

```@example landau
alpha = 0.1
kx = 0.5
nx = 50
np = 10000 * nx
mesh = OneDGrid( 0, 4π, nx)
x = LinRange(0, 4pi, nx+1)[1:end-1]
rng = MersenneTwister(42)
pa = landau_damping(rng, mesh, np, alpha, kx )
ex = zeros(Float64, nx)
rho = zeros(Float64, nx)
poisson = OneDPoisson( mesh )
pm = ParticleMeshCoupling(pa, mesh)
mat = compute_coeffs(pm, pa)
compute_rho!(rho, mat, mesh, pa)
solve!(ex, poisson, rho)
```

```@example landau
p = plot(layout=2)
plot!(p[1], x, ex)
plot!(p[1], x, alpha/kx * sin.(kx * x))
plot!(p[2], x, rho)
plot!(p[2], x, alpha * cos.(kx * x))
```

## Full simulation


```@example landau

using ParticleInCell

function main(nt, dt)
    
    nx = 50
    np = 10000 * nx
    mesh = OneDGrid( 0, 4π, nx)
    poisson = OneDPoisson( mesh )
    rng = MersenneTwister(42)
    α = 0.1
    kx = 0.5
    pa = landau_damping(rng, mesh, np, α, kx )
    pm = ParticleMeshCoupling(pa, mesh)
    energy = Float64[]
    e = zeros(Float64, nx)
    ρ = zeros(Float64, nx)
    xmin = mesh.xmin
    xmax = mesh.xmax
    mat = compute_coeffs(pm, pa)
    compute_rho!(ρ, mat, mesh, pa)
    solve!(e, poisson, ρ)
    for it in 1:nt+1       
        update_positions!(pa, mesh, dt)
        mat = compute_coeffs(pm, pa)
        compute_rho!(ρ, mat, mesh, pa)
        solve!(e, poisson, ρ)
        update_velocities!(pa, e, mat, dt)
        push!(energy, 0.5 * sum(e.^2) * mesh.dx) 
    end
    energy
end
```

```@example landau
nt, dt = 1000, 0.1
results = main(nt, dt)
t = collect(0:nt) .* dt
plot( t, results, yaxis = :log )
```


