# Vlasov-Poisson 

Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

```@example vp2d2v
using Plots
using ParticleInCell

dt = 0.1
nsteps = 100
alpha = 0.1
kx = 0.5

nx = 32
ny = 32
xmin = 0.0
ymin = 0.0
xmax = 4π
ymax = 4pi	

n_particles = 100000
degree_smoother = 3

mesh = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)

particles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)

sampler = LandauDamping( alpha, kx )

sample!( particles, sampler)

particles.array[5,:]  .= (4π * 4pi) ./ n_particles;
```

```@example vp2d2v
histogram( particles.array[1,:], normalized=true)
```

```@example vp2d2v
histogram( particles.array[2,:], normalized=true)
```

```@example vp2d2v
histogram( particles.array[3,:], normalized=true)
xlims!(-6,6)
```

```@example vp2d2v
histogram( particles.array[4,:], normalized=true)
xlims!(-6,6)
```

```@example vp2d2v
poisson = Poisson2DPeriodic( mesh )
kernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)

ex = zeros(nx, ny)
ey = zeros(nx, ny)
rho_dofs = zeros(nx*ny)

for i_part = 1:particles.n_particles
    xi = particles.array[1, i_part]
    yi = particles.array[2, i_part]
    wi = particles.array[5, i_part]
    ParticleInCell.add_charge!(rho_dofs, kernel, xi, yi, wi)
end
rho = reshape(rho_dofs, nx, ny )
surface(rho)
```

```@example vp2d2v
solve!(ex, ey, poisson, rho)
surface(ex)
```

```@example vp2d2v

poisson = Poisson2DPeriodic( mesh )

kernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)

problem = PICPoisson2D(  poisson, kernel )

dt = 0.1
nsteps = 100
alpha = 0.1
kx = 0.5


sampler = LandauDamping( alpha, kx )

sample!( particles, sampler)

particles.array[2,:] .*= ( ymax - ymin)
particles.array[5,:]  .= 4π * 4π / n_particles;

propagator = SplittingOperator( problem, particles ) 

energy = Float64[]

for j=1:100

    ParticleInCell.operator_t!(propagator, 0.5dt)
    ParticleInCell.charge_deposition!(propagator)
    ParticleInCell.solve_fields!(propagator)
    ParticleInCell.operator_v!(propagator, dt)
    ParticleInCell.operator_t!(propagator, 0.5dt)

    push!(energy, compute_field_energy(problem, 1))
          
end

plot(log.(energy))
```
