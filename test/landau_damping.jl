@testset "Landau Damping" begin

    nx = 128  # nombre de pts suivant x
    ny = 16   # nombre de pts suivant y
    alpha = 0.1
    kx = 0.5
    ky = 0.0
    dimx = 2pi / kx
    dimy = 1

    mesh = Mesh(dimx, nx, dimy, ny)

    dx = mesh.dx
    dy = mesh.dy

    fdtd = FDTD(mesh)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny
        fdtd.ex[i, j] = alpha / kx * sin(kx * mesh.x[i])
    end

    nbpart = 100 * nx * ny

    particles = Particles(nbpart)

    landau_sampling!(particles, alpha, kx)

    for istep = 1:10

        istep > 1 && faraday!(fdtd, 0.5dt)

        interpol_eb!(particles, fdtd)

        push_v!(particles, dt)

        push_x!(particles, mesh, 0.5dt)

        compute_current!(fdtd, particles)

        push_x!(particles, mesh, 0.5dt) 

        faraday!(fdtd, 0.5dt)

        ampere_maxwell!(fdtd, dt)

        time = time + dt

    end

    @test true

end
