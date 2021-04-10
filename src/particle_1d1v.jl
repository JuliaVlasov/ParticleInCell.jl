export TwoStreamInstability

struct TwoStreamInstability
    grid :: OneDGrid
end


function sample!(pg::ParticleGroup{1,1}, d::TwoStreamInstability)

    np = pg.n_particles
    xmin = d.grid.xmin
    xmax = d.grid.xmax
    dx = d.grid.dx

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

    pg.array[:,1] .= xp
    pg.array[:,2] .= vp

end

