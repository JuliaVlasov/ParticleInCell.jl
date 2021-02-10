using ParticleInCell
using LinearAlgebra


@testset "FDTD solver" begin

    dimx, dimy = 2π, 2π
    nx, ny = 128, 128
    dt = 1e-4
    nstep = 10

    mesh = Mesh(dimx, nx, dimy, ny)

    fdtd = FDTD(mesh)

    ω = sqrt(2)

    ex = zeros(nx, ny)
    ey = zeros(nx, ny)
    bz = zeros(nx, ny)
    jx = zeros(nx, ny)
    jy = zeros(nx, ny)

    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    x = mesh.x[1:nx]
    y = transpose(mesh.y[1:ny])
    xc = ( x[1:nx-1] .+ x[2:nx] ) ./ 2
    yc = ( y[1:ny-1] .+ y[2:ny] ) ./ 2

    t = 0
    ex .= 0
    ey .= 0
    t = 0.5dt
    bz .= -cos.(x) .* cos.(y) .* cos(ω * t)
    fdtd.bz .= cos.(xc) .* cos.(yc) .* cos(ω * t)

    sol_ex = +cos.(x) .* sin.(y) ./ ω
    sol_ey = -sin.(x) .* cos.(y) ./ ω
    sol_bz = -cos.(x) .* cos.(y)

    for istep = 1:nstep # Loop over time

        ampere_maxwell!(ex, ey, fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(ex .- sol_ex .* sin(ω * t))) < 0.001
        @test maximum(abs.(ey .- sol_ey .* sin(ω * t))) < 0.001

        faraday!(bz, fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(bz .- sol_bz .* cos(ω * t))) < 0.001

    end # next time step

end
