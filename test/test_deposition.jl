using Random
import ParticleInCell.F90

@testset "deposition" begin


    dimx, dimy = 4π, 4π
    nx, ny = 128, 128
    mesh = TwoDGrid(dimx, nx, dimy, ny)
    fdtd = FDTD(mesh)
    rng = MersenneTwister(1234)
    nbpart = 1_000_000
    jx = zeros(nx + 1, ny + 1)
    jy = zeros(nx + 1, ny + 1)

    kernel = CloudInCell()

    p = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    randn!(rng, p.array)
    p.array[1:2, :] .+= 2π
    p.array[3:4, :] .= (dimx * dimy) / nbpart
    p.array[5, :] .= (dimx * dimy) / nbpart

    @test sum(view(p.array, 1, :)) / nbpart ≈ 2π atol = 1e-2
    @test sum(view(p.array, 2, :)) / nbpart ≈ 2π atol = 1e-2

    compute_current!(jx, jy, mesh, kernel, p)

    @test sum(jx) ≈ 2.5872575761191685
    @test sum(jy) ≈ 2.5872575761191685

    F90.compute_current!(jx, jy, mesh, kernel, p)

    @test sum(jx) ≈ 2.5872575761191685
    @test sum(jy) ≈ 2.5872575761191685

end
