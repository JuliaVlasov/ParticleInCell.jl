using Plots
using ProgressMeter
using Random
using Revise

using ParticleInCell

const dt = 0.005     # Time step
const nt = 10000     # Number of time steps
const L  = 20π        #  Domain size 
const nx = 320       # Number of grid cells
const np = nx * 20   # Number of particles


mesh = Mesh1D( 0, 20π, nx)
rng = MersenneTwister(42)
poisson = Poisson1D( mesh )
particles = tsi(rng, mesh, np )
pm = ParticleMeshCoupling1D(particles, mesh)

function main()

    mesh = Mesh1D( 0, 20π, nx)
    poisson = Poisson1D( mesh )
    rng = MersenneTwister(42)
    particles = tsi(rng, mesh, np )
    pm = ParticleMeshCoupling1D(particles, mesh)
    b = Progress(nt+1)
    energy = Float64[]
    e = zeros(Float64, nx)
    ρ = zeros(Float64, nx)
    
    for it in 1:nt+1
        update_positions!(particles, mesh, dt) 
        mat = compute_coeffs(pm, particles)
        compute_rho!(ρ, mat, mesh, particles)
        solve!(e, poisson, ρ)
        update_velocities!(particles, e, mat, dt)
        push!(energy, 0.5 * sum(e.^2) * mesh.dx) 
        next!(b)
    end
    energy
end

@time results = main();

results



# +
t = (0:nt) .* dt

plot( t, results, yaxis=:log)
# -


