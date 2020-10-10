# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.6.0
#   kernelspec:
#     display_name: Julia 1.0.1
#     language: julia
#     name: julia-1.0
# ---

using Plots
using LinearAlgebra
using SparseArrays
using Random

# +
L  = 20π       # Domain size
DT = 0.005     # Time step
NT = 10000     # Number of time steps
NG = 320       # Number of grid cells
N  = NG * 20   # Number of particles
WP = 1.        # Plasma frequency
QM = -1.       # Charge/mass ratio
V0 = 0.9       # Stream velocity
VT = 0.0000001 # Thermal speed

# perturbation 
XP1 = 1.0 
mode = 1

Q = WP^2 / (QM*N/L)  # rho0*L/N: charge carried by a single particle?
rho_back = -Q*N/L  # Background charge density?
dx = L / NG # Grid step

# Auxilliary vectors
p = [1:N; 1:N]  # Some indices up to N
Poisson = spdiagm(-1 => ones(Float64,NG-2),
                   0 => -2*ones(Float64,NG),
                   1 => ones(Float64,NG-2))

# +
# Cell center coordinates
xg = collect(range(0, stop=L-dx, length=NG)) .+ dx/2
xp = collect(range(0, stop=L-L/N, length=N))  # Particle positions

rng = MersenneTwister(1234)

vp = VT .* (1 .- VT.^2).^(-0.5) .* randn(rng, Float64, N) # Particle momentum, initially Maxwellian
pm = collect(0:N-1)
pm = 1 .- 2 * mod.(pm.+1, 2)
vp .+= pm * (V0 * (1 - V0^2)^(-0.5)) # Momentum + stream velocity

# Add electron perturbation to excite the desired mode
xp .+= XP1 * (L/N) * sin.(2 * π * xp / L * mode)
xp[map(x -> x <  0, xp)] .+= L
xp[map(x -> x >= L, xp)] .-= L
# -

Poisson = spdiagm(-1 => ones(Float64,NG-2),
                   0 => -2*ones(Float64,NG),
                   1 => ones(Float64,NG-2))

function compute_rho(xp)
    # Project particles -> grid
    g1 = floor.(Int64,floor.(xp/dx .- 0.5))
    g = floor.(Int64,vcat(g1, g1.+1))
    fraz1 = 1 .- abs.(xp/dx .- g1 .- 0.5)
    fraz = vcat(fraz1, 1 .- fraz1)
    g[map( x -> x < 0   , g)] .+= NG
    g[map( x -> x > NG-1, g)] .-= NG
    g .+= 1
    println(minimum(p)," - ", maximum(p))
    println(minimum(g)," - ", maximum(g))
    mat = sparse(p, g, fraz, N, NG)
    Q / dx * vec(sum(mat,dims=1)) .+ rho_back
end

function solve_poisson( rho )
    # Compute electric field potential
    Φ = zeros(Float64,NG)
    rho .*= (-dx^2)
    sol = Poisson \  rho
    Φ[1:end-1] .= sol[1:end-1]
    Φ[end] = 0.0
    # Electric field on the grid
    (circshift(Φ, 1) - circshift(Φ, -1)) ./ (2dx)
end  

rho = compute_rho(xp)  
@show typeof(rho)
Eg = solve_poisson(rho)
plot(xg, Eg)

# ![](Eg.png)
