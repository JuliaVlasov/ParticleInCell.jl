# This is a temporary implementation of PIC 1D simulation
# It will be remove and replaced
# I keep it just for benchmark purpose


export Particles1D

struct Particles1D

    np::Int
    qm::Float64
    xp::Vector{Float64}
    vp::Vector{Float64}
    qp::Float64
 
end

export tsi

function tsi(rng, mesh, np)

    xmin = mesh.xmin
    xmax = mesh.xmax
    dx = mesh.dx

    wp = +1.0  # plasma frequency
    qm = -1.0  # charge/mass ratio
    qp = wp^2 / (qm * np / (xmax - xmin))
    xp = collect(LinRange(xmin, xmax, np + 1))[1:end-1]  # Particle positions
    v0 = 0.9       # Stream velocity
    vt = 0.0000001 # Thermal speed

    # Particle momentum, initially Maxwellian
    vp = vt .* (1 .- vt .^ 2) .^ (-0.5) .* randn(rng, Float64, np)
    pm = collect(0:np-1)
    pm = 1 .- 2 * mod.(pm .+ 1, 2)
    vp .+= pm .* (v0 * (1 - v0^2)^(-0.5)) # Momentum + stream velocity

    # Add electron perturbation to excite the desired mode perturbation 
    xp1 = 1.0
    mode = 1
    xp .+= xp1 * dx * sin.(2π .* xp ./ (xmax - xmin) .* mode)
    xp .= mod.(xp, xmax - xmin)

    return Particles1D(np, qm, xp, vp, qp)

end

export landau_damping

function landau_damping(rng, mesh, np, alpha, kx)

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1 - x0) > 1e-12)
            p = x0 + alpha * sin(kx * x0) / kx
            f = 1 + alpha * cos(kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end

    xp = zeros(np)
    vp = zeros(np)
    s = Sobol.SobolSeq(2)

    for i = 1:np
        v = sqrt(-2 * log((i - 0.5) / np))
        r1, r2 = Sobol.next!(s)
        θ = r1 * 2π
        xp[i] = newton(r2)
        vp[i] = v * cos(θ)
    end

    qm = 1.0
    qp = (mesh.xmax - mesh.xmin) / np

    return Particles1D(np, qm, xp, vp, qp)

end

export OneDPoisson

struct OneDPoisson

    dx::Any
    matrix::Any

    function OneDPoisson(mesh)

        nx = mesh.nx

        matrix = spdiagm(
            -1 => ones(Float64, nx - 2),
            0 => -2 * ones(Float64, nx),
            1 => ones(Float64, nx - 2),
        )
        new(mesh.dx, matrix)

    end

end

export solve!

"""
Compute electric field from charge density
"""
function solve!(e, poisson, ρ)
    dx = poisson.dx
    ρ .*= (-dx^2)
    ρ .= poisson.matrix \ ρ
    e[1:end-1] .= ρ[1:end-1]
    e[end] = 0.0
    e .= (circshift(e, 1) .- circshift(e, -1)) ./ (2dx)
end

export ParticleMeshCoupling

struct ParticleMeshCoupling

    np::Any
    nx::Any
    dx::Any
    g::Vector{Int}
    f::Vector{Float64}
    p::Any

    function ParticleMeshCoupling(p, mesh)
        np = p.np
        g = zeros(Int, np)
        f = zeros(Float64, np)
        p = [1:np; 1:np]
        new(np, mesh.nx, mesh.dx, g, f, p)
    end

end

export compute_coeffs

function compute_coeffs(pm, p)
    dx = pm.dx
    pm.f .= p.xp ./ dx .- 0.5
    pm.g .= floor.(Int, pm.f)
    g = vcat(pm.g, pm.g .+ 1)
    pm.f .= 1 .- abs.(pm.f .- pm.g)
    f = vcat(pm.f, 1 .- pm.f)
    g .= mod1.(g, pm.nx)
    mat = sparse(pm.p, g, f, pm.np, pm.nx)
    dropzeros!(mat)
    return mat
end

export compute_rho!

function compute_rho!(ρ, coeffs, mesh, p::Particles1D)

    xmin = mesh.xmin
    xmax = mesh.xmax
    dx = mesh.dx
    ρ_back = -p.qp * p.np / (xmax - xmin)
    ρ .= p.qp ./ dx .* vec(sum(coeffs, dims = 1)) .+ ρ_back

end

export update_positions!

"""
update particle position xp
"""
function update_positions!(p, mesh, dt)

    xmin = mesh.xmin
    xmax = mesh.xmax

    p.xp .+= p.vp .* dt
    p.xp .= xmin .+ mod.(p.xp .- xmin, xmax - xmin)

end

export update_velocities!

"""
update particle velocities vp
"""
function update_velocities!(p, e, coeffs, dt)

    p.vp .+= coeffs * e .* p.qm .* dt

end
