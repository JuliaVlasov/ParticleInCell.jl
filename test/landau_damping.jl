@testset "Landau Damping Julia" begin

    nx = 128
    ny = 16
    alpha = 0.1
    kx = 0.5
    ky = 0.0
    dimx = 2pi / kx
    dimy = 1.

    mesh1 = Mesh(dimx, nx, dimy, ny)
    mesh2 = Mesh(dimx, nx, dimy, ny)

    fdtd1 = FDTD(mesh1)
    fdtd2 = FDTD(mesh2)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny
        fdtd1.ex[i, j] = alpha / kx * sin(kx * mesh1.x[i])
        fdtd2.ex[i, j] = alpha / kx * sin(kx * mesh2.x[i])
    end

    nbpart = 100 * nx * ny

    particles1 = zeros(7, nbpart)
    particles2 = zeros(7, nbpart)

    landau_sampling!(particles1, alpha, kx)
    landau_sampling!(particles2, alpha, kx)

    for istep = 1:1

        istep > 1 && faraday!(fdtd1, mesh1, 0.5dt)
        istep > 1 && faraday!(fdtd2, mesh2, 0.5dt)

        update_fields!(mesh1, fdtd1)
        update_fields!(mesh2, fdtd2)

        interpolation!(particles1, mesh1)
        f90_interpolation!(particles2, mesh2)
        @test all(particles1 .== particles2)

        push_v!(particles1, dt)
        f90_push_v!(particles2, dt)
        @test all(particles1 .== particles2)

        push_x!(particles1, mesh1, 0.5dt)
        f90_push_x!(particles2, mesh2, 0.5dt)
        @test all(particles1 .== particles2)

        compute_current!(mesh1, particles1)
        f90_compute_current!(mesh2, particles2)

        push_x!(particles1, mesh1, 0.5dt) 
        push_x!(particles2, mesh2, 0.5dt) 

        faraday!(fdtd1, mesh1, 0.5dt)
        faraday!(fdtd2, mesh2, 0.5dt)

        update_currents!(fdtd1, mesh1)
        update_currents!(fdtd2, mesh2)

        ampere_maxwell!(fdtd1, mesh1, dt)
        ampere_maxwell!(fdtd2, mesh2, dt)

        time = time + dt

    end

    @test true

end
