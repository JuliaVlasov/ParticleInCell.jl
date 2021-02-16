@testset "Landau Damping Julia" begin

    nx = 128
    ny = 16
    alpha = 0.1
    kx = 0.5
    ky = 0.0
    dimx = 2pi / kx
    dimy = 1.

    mesh = Mesh(dimx, nx, dimy, ny)

    dx = mesh.dx
    dy = mesh.dy

    fdtd1 = FDTD(mesh)
    fdtd2 = FDTD(mesh)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny
        fdtd1.ex[i, j] = alpha / kx * sin(kx * mesh.x[i])
        fdtd2.ex[i, j] = alpha / kx * sin(kx * mesh.x[i])
    end

    nbpart = 100 * nx * ny

    particles1 = zeros(7, nbpart)
    particles2 = zeros(7, nbpart)

    landau_sampling!(particles1, nbpart, alpha, kx)
    landau_sampling!(particles2, nbpart, alpha, kx)

    for istep = 1:1

        istep > 1 && faraday!(fdtd1, 0.5dt)
        istep > 1 && faraday!(fdtd2, 0.5dt)

        interpol_eb!(particles1, nbpart, fdtd1)
        f90_interpolation!(particles2, fdtd2)
        @test all(particles1 .== particles2)

        push_v!(particles1, nbpart, dt)
        f90_push_v!(particles2, nbpart, dt)
        @test all(particles1 .== particles2)

        push_x!(particles1, nbpart, mesh, 0.5dt)
        f90_push_x!(particles2, nbpart, dimx, dimy, 0.5dt)
        @test all(particles1 .== particles2)

        @show maximum(fdtd1.ebj .- fdtd2.ebj)
        f90_deposition!(fdtd2, particles2)
        compute_current!(fdtd1, particles1, nbpart)
        @show maximum(fdtd1.ebj[1,:,:] .- fdtd2.ebj[1,:,:])
        @show maximum(fdtd1.ebj[2,:,:] .- fdtd2.ebj[2,:,:])
        @show maximum(fdtd1.ebj[3,:,:] .- fdtd2.ebj[3,:,:])
        @show maximum(fdtd1.ebj[4,:,:] .- fdtd2.ebj[4,:,:])
        @show maximum(fdtd1.ebj[5,:,:] .- fdtd2.ebj[5,:,:])
        @show maximum(fdtd1.jx .- fdtd2.jx)
        @show maximum(fdtd1.jy .- fdtd2.jy)

        push_x!(particles1, nbpart, mesh, 0.5dt) 

        faraday!(fdtd1, 0.5dt)

        ampere_maxwell!(fdtd1, dt)

        time = time + dt

    end

    @test true

end
