using Random
import ParticleInCell.F90

@testset "push" begin

    dimx, dimy = 1.0, 1.0
    nx, ny = 100, 100
    mesh = TwoDGrid(dimx, nx, dimy, ny)
    e = 1.0
    b = 2.0
    nbpart = 1
    pj = zeros(7, 1)
    dt = 0.1
    x0, y0 = 0.1, 0.1

    pj = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)
    pf = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    pj.array[1, 1] = x0
    pj.array[2, 1] = y0
    pj.array[5, 1] = 1.0

    ex = ones(nx + 1, ny + 1)
    ey = zeros(nx + 1, ny + 1)
    bz = ones(nx + 1, ny + 1)

    pf.array .= copy(pj.array)

    time = 0

    kernel = CloudInCell()

    for istep = 1:50

        push_v!(pj, kernel, mesh, ex, ey, bz, dt)
        push_x!(pj, mesh, dt)

        F90.push_v!(pf, kernel, mesh, ex, ey, bz, dt)
        F90.push_x!(pf, mesh, dt)

        @test all(pj.array .== pf.array)
    end

end
