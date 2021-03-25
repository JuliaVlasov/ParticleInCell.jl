@testset "Poisson 2D" begin


    eta1_min = 0.0
    eta1_max = 2π
    eta2_min = 0.0
    eta2_max = 2π

    nc_eta1 = 128
    nc_eta2 = 128

    ex = zeros(Float64, nc_eta1, nc_eta2)
    ey = zeros(Float64, nc_eta1, nc_eta2)

    rho = zeros(Float64, nc_eta1, nc_eta2)
    phi = zeros(Float64, nc_eta1, nc_eta2)

    ex_exact = zeros(nc_eta1, nc_eta2)
    ey_exact = zeros(nc_eta1, nc_eta2)
    phi_exact = zeros(nc_eta1, nc_eta2)

    grid = TwoDGrid(eta1_min, eta1_max, nc_eta1, eta2_min, eta2_max, nc_eta2)

    poisson = TwoDPoissonPeriodic(grid)

    mode = 2
    for i = 1:nc_eta1, j = 1:nc_eta2
        x1 = eta1_min + (i - 1) * (eta1_max - eta1_min) / nc_eta1
        x2 = eta2_min + (j - 1) * (eta2_max - eta2_min) / nc_eta2
        phi_exact[i, j] = -mode * sin(mode * x1) * cos(mode * x2)
        ex_exact[i, j] = mode^2 * cos(mode * x1) * cos(mode * x2)
        ey_exact[i, j] = -mode^2 * sin(mode * x1) * sin(mode * x2)
        rho[i, j] = -2 * mode^3 * sin(mode * x1) * cos(mode * x2)
    end

    rhs = copy(rho)
    solve!(phi, poisson, rhs)
    @test phi_exact ≈ phi

    rhs = copy(rho)
    solve!(ex, ey, poisson, rhs)
    @test ex_exact ≈ ex
    @test ey_exact ≈ ey

    x1_min = 0
    x1_max = 1

    x2_min = 0
    x2_max = 1

    nc_x1 = 32
    nc_x2 = 64

    phi = zeros(nc_x1, nc_x2)
    e1 = zeros(nc_x1, nc_x2)
    e2 = zeros(nc_x1, nc_x2)
    rho = ones(nc_x1, nc_x2)

    grid = TwoDGrid(x1_min, x1_max, nc_x1, x2_min, x2_max, nc_x2)

    poisson = TwoDPoissonPeriodic(grid)

    solve!(phi, poisson, rho)

    solve!(e1, e2, poisson, rho)

    nx, ny = 32, 32
    xmin, xmax = 0.0, 2π
    ymin, ymax = 0.0, 2π
    degree_smoother = 3
    n_particles = 1000

    particles = ParticleGroup{2,2}(n_particles, charge = 1.0, mass = 1.0, n_weights = 1)

    grid = TwoDGrid(xmin, xmax, nx, ymin, ymax, ny)

    poisson = TwoDPoissonPeriodic(grid)

    # Initialize the kernel smoother
    kernel_smoother = ParticleMeshCoupling2D(particles, grid, degree_smoother, :collocation)

    @test true

    # Initialize the PIC field solver
    pic_poisson = (poisson, kernel_smoother)

    @test true

end
