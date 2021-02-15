using ParticleInCell

@testset "Particles" begin

    nx = 64 # nombre de pts suivant x
    ny = 16   # nombre de pts suivant y

    alpha = 0.1
    kx = 0.5
    dimx = 2 * pi / kx
    dimy = 1

    mesh = Mesh(dimx, nx, dimy, ny)

    dx = mesh.dx
    dy = mesh.dy

    ex = zeros(nx, ny)
    rho = zeros(nx, ny)

    for i = 1:nx
        aux1 = alpha / kx * sin(kx * mesh.x[i])
        aux2 = alpha * cos(kx * mesh.x[i])
        for j = 1:ny
            ex[i, j] = aux1
            rho[i, j] = aux2
        end
    end

    nbpart = 200 * nx * ny

    particles = zeros(7, nbpart)

    landau_sampling!(particles, nbpart, alpha, kx)


    @test maximum(abs.(rho .- compute_rho(particles, mesh))) < 1e-2

end
