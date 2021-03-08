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

using Plots, Statistics
using Random
using Pkg
Pkg.activate(@__DIR__)
using Revise

# +
using ParticleInCell

alpha = 0.1
kx = 0.5
nx = 50
np = 10000 * nx
xmin, xmax = 0, 4π
mesh = OneDGrid( xmin, xmax, nx)
particles = ParticleGroup{1,1}(
        np;
        charge = 1.0,
        mass = 1.0,
        n_weights = 1)

sampler = LandauDamping( alpha, kx )
sample!( particles, mesh, sampler)
# -

p = plot(layout=2)
histogram!(p[1], particles.array[1,:], normalize=true, bins=50)
histogram!(p[2], particles.array[2,:], normalize=true, bins=50)

# +
function compute_rho!( rho, pm, particles)
    
    fill!(rho, 0.0)
 
    for i_part = 1:particles.n_particles
       xi = particles.array[1, i_part]
       wi = particles.array[3, i_part]
       add_charge!(rho, pm, xi, wi)
    end
    
end
# -

rho = zeros(Float64, nx)
compute_rho!(rho, pm, particles)

# +
maxwell = Maxwell1DFEM(mesh, spline_degree);

ex = zeros(Float64, nx)
x = LinRange(xmin, xmax, nx+1)[1:end-1]
compute_e_from_rho!( ex, maxwell, rho ) 
sval = eval_uniform_periodic_spline_curve(spline_degree-1, ex);
# -

p = plot(layout=2)
plot!(p[1], x, rho)
plot!(p[1], x, 1 .+ alpha * cos.(kx * x))
plot!(p[2], x, sval)
plot!(p[2], x, alpha/kx * sin.(kx * x))





