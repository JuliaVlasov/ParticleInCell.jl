@testset "CIC" begin

    nx = 64 # nombre de pts suivant x
    ny = 16   # nombre de pts suivant y

    alpha = 0.1
    kx = 0.5
    dimx = 2 * pi / kx
    dimy = 1

    mesh = TwoDGrid(dimx, nx, dimy, ny)

    dx = mesh.dx
    dy = mesh.dy

    rho = zeros(nx, ny)

    for i = 1:nx
        aux2 = alpha * cos(kx * mesh.x[i])
        for j = 1:ny
            rho[i, j] = aux2
        end
    end

    nbpart = 1000 * nx * ny

    p = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    sampler = LandauDamping(alpha, kx)

    sample!(p, mesh, sampler)

    kernel = CloudInCell()

    ρ = compute_rho(p, kernel, mesh)

    @test ρ ≈ rho atol=1e-2

end


