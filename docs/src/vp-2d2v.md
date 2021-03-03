# Vlasov-Poisson 

Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

```julia
using Pkg, Plots
Pkg.activate("/Users/navaro/JuliaProjects/ParticleInCell.jl/")
using Revise
```

```julia
using ParticleInCell

delta_t = 0.1
n_time_steps = 100
alpha = 0.1
kx = 1.0

ng_x1 = 32
ng_x2 = 32
x1_min = 0.0
x2_min = 0.0
x1_max = 4π
x2_max = 4π

n_particles = 100000
degree_smoother = 1

mesh = TwoDGrid( x1_min, x1_max, ng_x1, x2_min, x2_max, ng_x2)

particles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)

sampler = LandauDamping( alpha, kx )

sample!( particles, sampler)
```

```julia
histogram( particles.array[1,:], normalized=true)
```

```julia
histogram( particles.array[2,:], normalized=true)
```

```julia
histogram( particles.array[3,:], normalized=true)
xlims!(-6,6)
```

```julia
histogram( particles.array[4,:], normalized=true)
xlims!(-6,6)
```

```julia
poisson = Poisson2DPeriodic( mesh )
kernel = ParticleMeshCoupling2D( particles, mesh, 1, :collocation)

ex = zeros(ng_x1, ng_x2)
ey = zeros(ng_x1, ng_x2)
rho_dofs = zeros(ng_x1*ng_x2)

for i_part = 1:particles.n_particles
    xi = particles.array[1, i_part]
    yi = particles.array[2, i_part]
    wi = particles.array[5, i_part]
    ParticleInCell.add_charge!(rho_dofs, kernel, (xi, yi), wi)
end
```

```julia
solve!(ex, ey, poisson, rho)
```

```julia
kernel = ParticleMeshCoupling( particles, mesh, degree_smoother, :collocation)

problem = PicPoisson2DPeriodic( kernel, poisson )

propagator = SplittingOperator( pic ) 

energy = Float64[]

for j=1:nsteps

    strang_splitting!( propagator, delta_t)

    push!(energy, compute_field_energy(pic))
          
end

energy
```
