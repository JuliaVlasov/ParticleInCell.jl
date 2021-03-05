using ProgressMeter
using TimerOutputs
using ParticleInCell

const to = TimerOutput()

function run( nsteps)

    dt = 0.1
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
    
    sample!( particles, mesh, sampler)
    
    poisson = Poisson2DPeriodic( mesh )
    kernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)
    
    problem = PICPoisson2D(  poisson, kernel )
    
    particles.array[2,:] .*= ( ymax - ymin)
    particles.array[5,:]  .= (xmax - xmin) * (ymax - ymin) / n_particles;
    
    propagator = SplittingOperator( problem, particles ) 
    
    time = Float64[0.0]
    energy = Float64[compute_field_energy(problem, 1)]
    
    reset_timer!()
    @showprogress 1 for j=1:nsteps
    
        @timeit to "operatorT" ParticleInCell.operator_t!(propagator, 0.5dt)
        @timeit to "charge computation" ParticleInCell.charge_deposition!(propagator)
        @timeit to "fields solver" ParticleInCell.solve_fields!(propagator)
        @timeit to "operatorV" ParticleInCell.operator_v!(propagator, dt)
        @timeit to "operatorT" ParticleInCell.operator_t!(propagator, 0.5dt)
    
        push!(time, j*dt)
        push!(energy, compute_field_energy(problem, 1))
              
    end

    time, energy

end

t, energy = run( 1 ) # trigger building
@show nstep = 1000
@time t, energy = run( nstep )
show(to)

open("results_jl.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
