@testset "Operator Splitting" begin

    # Set parameters
    n_particles = 2
    eta_min = 0.0
    eta_max = 4π
    num_cells = 10
    delta_t = 0.1
    degree_smoother = 3

    pg = ParticleGroup{2,2}(n_particles, charge = 1.0, mass = 1.0, n_weights = 1)


    particle_info_ref = [
        11.780972450961723 5.4977871437821380
        0.78539816339744828 7.0685834705770345
        0.15731068461017067 -1.5341205443525459
        1.5341205443525459 -0.15731068461017067
        86.251495608834688 71.662174808595040
    ]

    for i_part = 1:n_particles
        pg.array[1, i_part] = particle_info_ref[1, i_part]
        pg.array[2, i_part] = particle_info_ref[2, i_part]
        pg.array[3, i_part] = particle_info_ref[3, i_part]
        pg.array[4, i_part] = particle_info_ref[4, i_part]
        pg.array[5, i_part] = particle_info_ref[5, i_part]
    end

    @test pg.array ≈ particle_info_ref

    grid = TwoDGrid(eta_min, eta_max, num_cells, eta_min, eta_max, num_cells)

    kernel = ParticleMeshCoupling2D(pg, grid, degree_smoother, :collocation)

    poisson = Poisson2DPeriodic(grid)

    pic = PICPoisson2D(poisson, kernel)

    propagator = SplittingOperator(pic, pg)

    ParticleInCell.operator_t!(propagator, delta_t)

    particle_info_ref = [
        11.796703519422740 5.3443750893468831
        0.93881021783270291 7.0528524021160175
        0.15731068461017067 -1.5341205443525459
        1.5341205443525459 -0.15731068461017067
        86.251495608834688 71.662174808595040
    ]

    @test pg.array ≈ particle_info_ref

    charge_deposition!(propagator)
    solve_fields!(propagator)
    ParticleInCell.operator_v!(propagator, delta_t)


    particle_info_ref = [
        11.796703519422740 5.3443750893468831
        0.93881021783270291 7.0528524021160175
        0.15346588344558551 -1.5294930005799796
        1.5302579166423358 -0.15266168508042444
        86.251495608834688 71.662174808595040
    ]

    @test pg.array ≈ particle_info_ref

end
