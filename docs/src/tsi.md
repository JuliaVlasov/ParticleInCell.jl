# Two-stream instability

```@setup tsi
using Plots
using Random
```

```@example tsi
using ParticleInCell

const dt = 0.005     # Time step
const nt = 10000     # Number of time steps
const L  = 20π       #  Domain size 
const nx = 320       # Number of grid cells
const np = nx * 20   # Number of particles


mesh = OneDGrid( 0, 20π, nx)
rng = MersenneTwister(42)
poisson = OneDPoisson( mesh )
particles = tsi(rng, mesh, np )
pm = ParticleMeshCoupling(particles, mesh)
```

```@example tsi
function main()

    mesh = OneDGrid( 0, 20π, nx)
    poisson = OneDPoisson( mesh )
    rng = MersenneTwister(42)
    pa = tsi(rng, mesh, np )
    pm = ParticleMeshCoupling(pa, mesh)
    energy = Float64[]
    e = zeros(Float64, nx)
    ρ = zeros(Float64, nx)
    xmin = mesh.xmin
    xmax = mesh.xmax
    
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

```@example tsi
results = main()
t = (0:nt) .* dt
plot( t, results, yaxis=:log)
```
