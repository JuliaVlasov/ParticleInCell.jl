using FFTW

export OneDPoisson

"""
    OneDPoisson

Derived type to solve the Poisson equation on 1d regular 
cartesian mesh with periodic boundary conditions
- kx   : wave number in x
- nc_x : cells number in x
- dx   : x step size
- rht  : fft(rho)
"""
struct OneDPoisson

    grid::OneDGrid
    kx::Array{Float64,1}
    rht::Array{ComplexF64,1}

    function OneDPoisson(grid::OneDGrid)

        nx = grid.nx
        rht = zeros(ComplexF64, div(nx, 2) + 1)
        kx = zeros(nx ÷ 2 + 1)

        kx0 = 2π / grid.dimx

        for ik = 1:nx÷2+1
            kx[ik] = (ik - 1) * kx0
        end

        kx[1] = 1.0

        new(grid, kx, rht)

    end

end


export compute_e_from_rho!

"""
    compute_e_from_rho!( ex, poisson, rho )

computes electric field `e` from `rho` by solving Poisson equation.

"""
function compute_e_from_rho!(ex, poisson::OneDPoisson, rho)

    poisson.rht .= rfft(rho)

    poisson.rht[1] = 0.0
    poisson.rht .*= -1im .* poisson.kx

    ex .= irfft(poisson.rht, poisson.grid.nx)

end

export OneDPoissonPIC

"""
    OneDPoissonPIC( poisson, kernel )

- kernel : Kernel smoother taking care of charge deposition and field evaluation
- poisson : Poisson solver
- rho : Coefficients of expansion of rho (MPI global version)
- efield : Coefficients of expansion of electric field
- phi : Coefficients of expansion of potential
"""
struct OneDPoissonPIC

    ndofs::Int
    kernel::ParticleMeshCoupling2D
    poisson::OneDPoisson
    rho::Vector{Float64}
    efield::Vector{Vector{Float64}}
    phi::Vector{Float64}

    function OneDPoissonPIC(poisson::OneDPoisson, kernel::ParticleMeshCoupling1D)

        ndofs = poisson.grid.nx

        rho = zeros(ndofs)
        efield = [zeros(ndofs), zeros(ndofs)]

        new(ndofs, kernel, poisson, rho, efield)

    end

end

"""
    add_charge!(pic, position, marker_charge)

Add charge from one particle
- self : Pic Poisson solver object
- position : Position of the particle
- marker_charge : Particle weight times charge
"""
function add_charge!(pic::OneDPoissonPIC, position::Float64, marker_charge)

    GEMPIC.add_charge!(pic.rho, pic.kernel, position, marker_charge)

end

"""
    evaluate_rho!(pic, position)

Evaluate charge density at rho at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of rho at given position
"""
function evaluate_rho!(pic :: OneDPoissonPIC, position)

    GEMPIC.evaluate!(kernel, position, pic.rho)

end

export OneDSplittingOperator

"""
Operator splitting type for 2d2v Vlasov-Poisson
- pic :: PIC poisson solver
- pg :: Particle group
"""
struct OneDSplittingOperator

    pic::OneDPoissonPIC
    pg::ParticleGroup{1,1}

end

"""
    operator_t(split, dt)
Push x 
- split :: time splitting object 
- dt   :: time step
"""
function operator_t!(split::OneDSplittingOperator, dt)

    for i_part = 1:split.pg.n_particles
        split.pg.array[1, i_part] += dt * split.pg.array[2, i_part]
    end

end

function operator_v!(split::OneDSplittingOperator, dt)

    qm = split.pg.q_over_m

    for i_part = 1:split.pg.n_particles

        # Evaluate efields at particle position

        xi = split.pg.array[1, i_part]

        ex = GEMPIC.evaluate(split.pic.kernel, xi, split.pic.efield)

        vx_new = split.pg.array[1, i_part]

        vx_new += dt * qm * ex

        split.pg.array[2, i_part] = vx_new

    end

end

export charge_deposition!

function charge_deposition!(split::OneDSplittingOperator)

    fill!(split.pic.rho_dofs, 0.0)
    for i_part = 1:split.pg.n_particles
        xp = split.pg.array[1, i_part]
        wp = split.pg.array[3, i_part]
        GEMPIC.add_charge!(split.pic.rho_dofs, split.pic.kernel, xp, wp)
    end

end

export solve_fields!

"""
    solve_fields!( pic )

Solve efields from rho
"""
function solve_fields!(pic :: OneDPoissonPIC)

    compute_e_from_rho!(pic.efield[1], pic.poisson, pic.rho)

end


"""
Solve Poisson's equation for the electric field
"""
solve_fields!(split::OneDSplittingOperator) = solve_fields!(split.pic)


export strang_splitting!

"""
    strang_splitting!(split, dt)

Strang splitting
- split :: time splitting object 
- dt   :: time step
"""
function strang_splitting!(split::OneDSplittingOperator, dt)

    operator_t!(split, 0.5dt)
    charge_deposition!(split)
    solve_fields!(split)
    operator_v!(split, dt)
    operator_t!(split, 0.5dt)

end

export compute_field_energy

"""
    compute_field_energy

Compute the squared l2 norm of electric field
"""
compute_field_energy(pic::OneDPoissonPIC, component::Int) =
    sum(pic.efield[component] .^ 2) * pic.poisson.grid.dx
