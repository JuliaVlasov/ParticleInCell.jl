# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

using ParticleInCell
using Plots

# +
nx = 64 # nombre de pts suivant x
ny = 16   # nombre de pts suivant y

alpha = 0.1
kx = 0.5
dimx = 2 * pi / kx
dimy = 1

mesh = TwoDGrid(dimx, nx, dimy, ny)

dx = mesh.dx
dy = mesh.dy

rho = zeros(nx, ny)

for i = 1:nx
    aux2 = alpha * cos(kx * mesh.x[i])
    for j = 1:ny
        rho[i, j] = aux2
    end
end
# -

surface(rho)

# +
nbpart = 10 * nx * ny

p = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

sampler = LandauDamping(alpha, kx)

ParticleInCell.sample!(p, mesh, sampler)

kernel = CloudInCell()

ρ = compute_rho(p, kernel, mesh)
# -

surface(ρ)

# +
poisson = TwoDPoissonPeriodic(mesh)

ϕ = similar(ρ)
solve!(ϕ, poisson, ρ)
# -

surface(ϕ)

# +
ex = similar(ρ)
ey = similar(ρ)

solve!(ex, ey, poisson, ρ)
# -

surface(ex)

surface(ey)

# +
function interpolate_cic!(epx, epy, mesh, p, ex, ey)

    dx = mesh.dx
    dy = mesh.dy

    for m = eachindex(epx, epy)
    
       xp = p.array[1, m] / dx
       yp = p.array[2, m] / dy

       i = trunc(Int, xp) + 1
       j = trunc(Int, yp) + 1

       dxp = xp - i + 1
       dyp = yp - j + 1

       a1 = (1 - dxp) * (1 - dyp)
       a2 = dxp * (1 - dyp)
       a3 = dxp * dyp
       a4 = (1 - dxp) * dyp

       ip1 = mod1(i + 1, nx)
       jp1 = mod1(j + 1, ny)

       epx[m] = a1 * ex[i,  j] + a2 * ex[ip1, j] + a3 * ex[ip1,  jp1] + a4 * ex[i, jp1]
       epy[m] = a1 * ey[i,  j] + a2 * ey[ip1, j] + a3 * ey[ip1,  jp1] + a4 * ey[i, jp1]
                
    
    end

end

# +
epx = zeros(p.n_particles)
epy = zeros(p.n_particles)

interpolate_cic!(epx, epy, mesh, p, ex, ey)
# -

scatter(p.array[1,:], p.array[2,:], epx, ms=1)

scatter(p.array[1,:], p.array[2,:], epy, ms=1)


