# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:nomarker
#     text_representation:
#       extension: .jl
#       format_name: nomarker
#       format_version: '1.0'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Julia 1.5.2
#     language: julia
#     name: julia-1.5
# ---

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
    pa = tsi(rng, mesh, np )
    pm = ParticleMeshCoupling1D(pa, mesh)
    b = Progress(nt+1)
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
        next!(b)
    end
    energy
end

@time results = main()
t = (0:nt) .* dt
plot( t, results, yaxis=:log)




