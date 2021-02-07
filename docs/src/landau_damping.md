# Landau damping


```@setup landau
using Plots
using Random
```


```@example landau

using ParticleInCell

function main(nt, dt)
    
    nx = 50
    np = 10000 * nx
    mesh = Mesh1D( 0, 4π, nx)
    poisson = Poisson1D( mesh )
    rng = MersenneTwister(42)
    α = 0.5
    kx = 0.5
    pa = landau_damping(rng, mesh, np, α, kx )
    pm = ParticleMeshCoupling1D(pa, mesh)
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
nt, dt = 2000, 0.01
results = main(nt, dt)
t = collect(0:nt) .* dt
plot( t, results, yaxis = :log )
```


