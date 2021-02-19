using ParticleInCell
using LinearAlgebra


@testset "FDTD solver" begin

    dimx, dimy = 2π, 2π
    nx, ny = 128, 128
    dt = 1e-4
    nstep = 8

    mesh = TwoDGrid(dimx, nx, dimy, ny)

    fdtd = FDTD(mesh)

    ω = sqrt(2)

    xn = mesh.x
    yn = mesh.y |> transpose
    xc = ( mesh.x[1:end-1] .+ mesh.x[2:end] ) ./ 2
    yc = ( mesh.y[1:end-1] .+ mesh.y[2:end] ) ./ 2 |> transpose

    # Ex and Ey are set to zero at t = 0.0
    # Bz is set at  t = dt/2
    t = 0.5dt
    fdtd.bz .= - cos.(xc) .* cos.(yc) .* cos(ω * t)

    sol_ex(x, y) = +cos.(x) .* sin.(y) ./ ω
    sol_ey(x, y) = -sin.(x) .* cos.(y) ./ ω
    sol_bz(x, y) = -cos.(x) .* cos.(y)

    for istep = 1:nstep # Loop over time

        ampere_maxwell!(fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(fdtd.ex .- sol_ex(xc,yn) .* sin(ω * t))) < 1e-6
        @test maximum(abs.(fdtd.ey .- sol_ey(xn,yc) .* sin(ω * t))) < 1e-6

        update_fields!(mesh, fdtd)

        @test maximum(abs.(mesh.ex .- sol_ex(xn,yn) .* sin(ω * t))) < 1e-6
        @test maximum(abs.(mesh.ey .- sol_ey(xn,yn) .* sin(ω * t))) < 1e-6

        faraday!(fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(fdtd.bz .- sol_bz(xc,yc) .* cos(ω * t))) < 1e-6

        update_fields!(mesh, fdtd)

        @test maximum(abs.(mesh.bz .- sol_bz(xn,yn) .* cos(ω * t))) < 1e-3

    end # next time step

end
