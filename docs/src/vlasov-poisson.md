# Vlasov-Poisson 

Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.

```@example vp2d2v
using Plots
using GEMPIC
using ParticleInCell

dt = 0.1
nsteps = 100
alpha = 0.1
kx = 0.5

nx = 128
ny = 16
xmin, xmax = 0.0, 4π
ymin, ymax = 0.0, 1.0

n_particles = 100000
degree_smoother = 3

mesh = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)

particles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)

sampler = LandauDamping( alpha, kx )

ParticleInCell.sample!( particles, mesh, sampler)

particles.array[5,:]  .= (xmax - xmin) * (ymax - ymin) ./ n_particles;
```

```@example vp2d2v
p = plot(layout=2)
histogram!( p[1], particles.array[1,:], normalized=true)
histogram!( p[2], particles.array[2,:], normalized=true)
```

```@example vp2d2v
p = plot(layout=2)

histogram!( p[1], particles.array[3,:], normalized=true)
xlims!(-6,6)
histogram!( p[2], particles.array[4,:], normalized=true)
xlims!(-6,6)
```


```@example vp2d2v
poisson = TwoDPoissonPeriodic( mesh )
kernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)

ex = zeros(nx, ny)
ey = zeros(nx, ny)
rho_dofs = zeros(nx*ny)

for i_part = 1:particles.n_particles
    xi = particles.array[1, i_part]
    yi = particles.array[2, i_part]
    wi = particles.array[5, i_part]
    GEMPIC.add_charge!(rho_dofs, kernel, xi, yi, wi)
end
rho = reshape(rho_dofs, nx, ny )
solve!(ex, ey, poisson, rho)
p = plot(layout=(2))
surface!(p[1], ex)
surface!(p[2],rho)
```

```@example vp2d2v

poisson = TwoDPoissonPeriodic( mesh )

kernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)

problem = TwoDPoissonPIC(  poisson, kernel )

dt = 0.1
nsteps = 100
alpha = 0.1
kx = 0.5

particles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)

sampler = LandauDamping( alpha, kx )

ParticleInCell.sample!( particles, mesh, sampler)

particles.array[2,:] .*= ( ymax - ymin)
particles.array[5,:]  .= (xmax - xmin) * (ymax - ymin) / n_particles;

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
