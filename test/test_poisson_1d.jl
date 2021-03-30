@testset "Poisson 1D" begin

    alpha = 0.1
    kx = 0.5
    xmin, xmax = 0.0, 2π / kx
    nx = 64

    mesh = OneDGrid( xmin, xmax, nx )

    rho = zeros(nx)
    expected_ex = zeros(nx)

    for i in 1:nx
        rho[i] = alpha * cos(kx * mesh.x[i])
        expected_ex[i] = alpha / kx * sin(kx * mesh.x[i])
    end

    poisson = OneDPoisson( mesh )
    
    computed_ex = zero(expected_ex)

    ParticleInCell.compute_e_from_rho!( computed_ex, poisson, rho )

    @show maximum(abs.( computed_ex .- expected_ex ))

    @test computed_ex ≈ expected_ex


end
