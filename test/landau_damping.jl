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

    maxwell = Maxwell(mesh)

    ex = zeros(nx, ny)
    ey = zeros(nx, ny)
    bz = zeros(nx, ny)
    jx = zeros(nx, ny)
    jy = zeros(nx, ny)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny
        ex[i, j] = alpha / kx * sin(kx * mesh.x[i])
    end

    nbpart = 100 * nx * ny

    particles = Particles(nbpart)

    landau_sampling!(particles, alpha, kx)
    update_cells!(particles, mesh)

    for istep = 1:10
        if istep > 1
            faraday!(bz, maxwell, ex, ey, 0.5dt)
        end

        interpol_eb!(ex, ey, bz, particles, mesh)
        push_v!(particles, dt)

        push_x!(particles, mesh, 0.5dt)
        compute_current!(jx, jy, particles, mesh)
        push_x!(particles, mesh, 0.5dt)  # x(n+1/2) -- x(n+1)
        faraday!(bz, maxwell, ex, ey, 0.5dt)
        ampere_maxwell!(ex, ey, maxwell, bz, jx, jy, dt)

        time = time + dt

    end

    @test true

end
