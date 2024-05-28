# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# # Landau damping

using ParticleInCell
using Plots
using ProgressMeter
using DispersionRelations

# +
nx = 128
ny = 16
alpha = 0.1
kx = 0.5
ky = 0.0
dimx = 2pi / kx
dimy = 1.0
dt = 0.01

mesh = TwoDGrid(dimx, nx, dimy, ny)
ρ = zeros(nx, ny)
for i = 1:nx, j = 1:ny
    ρ[i, j] = alpha * cos(kx * mesh.x[i])
end

surface(ρ)

# +
nbpart = 100 * nx * ny

p = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

sampler = LandauDamping(alpha, kx)

ParticleInCell.sample!(p, mesh, sampler)

kernel = CloudInCell()

ρ .= compute_rho(p, kernel, mesh)

surface(ρ)

# +
ex = similar(ρ)
ey = similar(ρ)

poisson = TwoDPoissonPeriodic(mesh)

solve!(ex, ey, poisson, ρ)

nsteps = 1000

# +
energy(ex) = sqrt(sum( ex .^ 2) * mesh.dx * mesh.dy)

t = [0.0]
e = [energy(ex)]

@showprogress 1 for i in 1:nsteps
    
    push_x!(p, mesh, 0.5dt)
    
    ρ = compute_rho(p, kernel, mesh)
    solve!(ex, ey, poisson, ρ)
    push!(t, i*dt)
    push!(e, energy(ex))
    
    push_v!(p, kernel, mesh, ex, ey, bz, dt)

    push_x!(p, mesh, 0.5dt)
    
end
# -
line, γ = fit_complex_frequency(t, e)
plot(t, e, yscale = :log10)  
plot!(t, line, label = "$(imag(γ))")