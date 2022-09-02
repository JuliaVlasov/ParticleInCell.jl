# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.7.1
#   kernelspec:
#     display_name: Julia 1.5.3
#     language: julia
#     name: julia-1.5
# ---

using Plots, Statistics, GEMPIC
using Random
using Pkg
Pkg.activate(@__DIR__)
using Revise

# +
using ParticleInCell

alpha = 0.1
kx = 0.5
nx = 64
np = 10000 * nx
xmin, xmax = 0, 2π/kx
mesh = OneDGrid( xmin, xmax, nx)
particles = ParticleGroup{1,1}(
        np;
        charge = 1.0,
        mass = 1.0,
        n_weights = 1)

sampler = LandauDamping( alpha, kx )
ParticleInCell.sample!( particles, mesh, sampler)

particles.array[3,:] .= (xmax - xmin) / np

# -

p = plot(layout=2)
histogram!(p[1], particles.array[1,:], normalize=true, bins=50)
histogram!(p[2], particles.array[2,:], normalize=true, bins=50)
display(p)

degree_smoother = 3

kernel = ParticleMeshCoupling1D( mesh, np, degree_smoother, :collocation)

# +
function compute_rho!( rho, kernel, particles)
    
    fill!(rho, 0.0)
 
    for i_part = 1:particles.n_particles
       xi = particles.array[1, i_part]
       wi = particles.array[3, i_part]
       GEMPIC.add_charge!(rho, kernel, xi, wi)
    end
    
end
# -

rho = zeros(Float64, nx)
compute_rho!(rho, kernel, particles)

# +
spline_degree = 3
maxwell = Maxwell1DFEM(mesh, spline_degree);

ex = zeros(Float64, nx)
x = LinRange(xmin, xmax, nx+1)[1:end-1]
compute_e_from_rho!( ex, maxwell, rho ) 
# -

p = plot(layout=2)
sval = eval_uniform_periodic_spline_curve(spline_degree-1, rho);
plot!(p[1], x, sval, label="computed")
plot!(p[1], x, 1 .+ alpha * cos.(kx * x), label="expected")
sval = eval_uniform_periodic_spline_curve(spline_degree-1, ex);
plot!(p[2], x, sval, label="computed")
plot!(p[2], x, alpha/kx * sin.(kx * x), label="expected")


