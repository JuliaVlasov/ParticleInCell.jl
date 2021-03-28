@testset "Poisson 1D" begin

    alpha = 0.1
    kx = 0.5
    xmin, xmax = 0.0, 2π / kx
    nx = 64

    grid = OneDGrid( xmin, xmax, nx )

    for i in 1:nx
        rho[i] = alpha * cos.(kx * x)
        ex[i] = alpha / kx * sin(kx * grid.x[i])
    end

    poisson = OneDPoisson( grid )
    
    @test true


end
