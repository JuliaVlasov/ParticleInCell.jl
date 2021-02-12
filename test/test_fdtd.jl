using ParticleInCell
using LinearAlgebra


@testset "FDTD solver" begin

    dimx, dimy = 2π, 2π
    nx, ny = 128, 128
    dt = 1e-4
    nstep = 8

    mesh = Mesh(dimx, nx, dimy, ny)

    fdtd = FDTD(mesh)

    ω = sqrt(2)

    eb = zeros(3, nx, ny)
    jx = zeros(nx, ny)
    jy = zeros(nx, ny)

    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    x = mesh.x[1:nx]
    y = mesh.y[1:ny] |> transpose
    xc = ( mesh.x[1:nx] .+ mesh.x[2:nx+1] ) ./ 2
    yc = ( mesh.y[1:ny] .+ mesh.y[2:ny+1] ) ./ 2 |> transpose

    t = 0
    fill!(eb, 0)
    t = 0.5dt
    fdtd.bz .= - cos.(xc) .* cos.(yc) .* cos(ω * t)

    sol_ex(x, y) = +cos.(x) .* sin.(y) ./ ω
    sol_ey(x, y) = -sin.(x) .* cos.(y) ./ ω
    sol_bz(x, y) = -cos.(x) .* cos.(y)

    for istep = 1:nstep # Loop over time

        ampere_maxwell!(eb, fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(fdtd.ex .- sol_ex(xc,y) .* sin(ω * t))) < 1e-6
        @test maximum(abs.(fdtd.ey .- sol_ey(x,yc) .* sin(ω * t))) < 1e-6
        @test maximum(abs.(eb[1,:,:] .- sol_ex(x,y) .* sin(ω * t))) < 1e-6
        @test maximum(abs.(eb[2,:,:] .- sol_ey(x,y) .* sin(ω * t))) < 1e-6

        faraday!(eb, fdtd, mesh, dt)

        t = t + 0.5dt

        @test maximum(abs.(fdtd.bz .- sol_bz(xc,yc) .* cos(ω * t))) < 1e-6
        @test maximum(abs.(eb[3,:,:] .- sol_bz(x,y) .* cos(ω * t))) < 1e-3

    end # next time step

end
