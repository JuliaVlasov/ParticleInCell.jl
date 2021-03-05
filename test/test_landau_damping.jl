import ParticleInCell.F90

@testset "Landau Damping Julia" begin

    nx = 128
    ny = 16
    alpha = 0.1
    kx = 0.5
    ky = 0.0
    dimx = 2pi / kx
    dimy = 1.0

    mesh1 = TwoDGrid(dimx, nx, dimy, ny)
    mesh2 = TwoDGrid(dimx, nx, dimy, ny)

    nbpart = 100 * nx * ny

    group1 = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)
    group2 = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    alpha = 0.1
    kx = 0.5

    sampler = LandauDamping(alpha, kx)

    sample!(group1, mesh1, sampler)
    sample!(group2, mesh2, sampler)

    fdtd1 = FDTD(mesh1)
    fdtd2 = FDTD(mesh2)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny+1
        fdtd1.ex[i, j] = alpha / kx * sin(kx * mesh1.x[i])
        fdtd2.ex[i, j] = alpha / kx * sin(kx * mesh2.x[i])
    end

    kernel = CloudInCell()

    particles1 = group1.array
    particles2 = group2.array

    for istep = 1:1

        istep > 1 && faraday!(fdtd1, mesh1, 0.5dt)
        istep > 1 && faraday!(fdtd2, mesh2, 0.5dt)

        update_fields!(mesh1, fdtd1)
        update_fields!(mesh2, fdtd2)

        push_v!(group1, kernel, mesh1,  dt)
        F90.push_v!(group2, kernel, mesh2, dt)
        @test all(particles1 .== particles2)

        push_x!(group1, mesh1, 0.5dt)
        F90.push_x!(group2, mesh2, 0.5dt)
        @test all(particles1 .== particles2)

        compute_current!(mesh1, kernel, group1)
        F90.compute_current!(mesh2, kernel, group2)

        push_x!(group1, mesh1, 0.5dt)
        F90.push_x!(group2, mesh2, 0.5dt)
        @test all(particles1 .== particles2)

        faraday!(fdtd1, mesh1, 0.5dt)
        faraday!(fdtd2, mesh2, 0.5dt)

        ampere_maxwell!(fdtd1, mesh1, dt)
        ampere_maxwell!(fdtd2, mesh2, dt)

        time = time + dt

    end

    @test true

end
