using Random
import ParticleInCell.F90

@testset "push"  begin

    dimx, dimy = 1.0, 1.0
    nx, ny = 100, 100
    mesh = TwoDGrid( dimx, nx, dimy, ny)
    e = 1.0
    b = 2.0
    nbpart = 1
    pj = zeros(7, 1)
    dt   = 0.1
    x0, y0 = 0.1, 0.1
    pj[1,1] = x0
    pj[2,1] = y0
    pj[5,1] = 1.0
    pj[7,1] = 1.0

    pf = copy(pj)

    time = 0

    for istep = 1:50
        push_v!( pj, dt)
        push_x!( pj, mesh, dt)

        F90.push_v!( pf, dt)
        F90.push_x!( pf, mesh, dt)

        @test all( pj .== pf)
    end

end
