using Random

@testset "push"  begin

    dimx, dimy = 1.0, 1.0
    nx, ny = 100, 100
    mesh = Mesh( dimx, nx, dimy, ny)
    e = 1.0
    b = 2.0
    nbpart = 1
    pj = zeros(7, 1)
    dt   = 1.0
    x0, y0 = 0.1, 0.1
    pj[1,1] = x0
    pj[2,1] = y0
    pj[5,1] = 1.0
    pj[7,1] = 1.0

    pf = copy(pj)

    time = 0

    for istep in 1:50
        
        push_v!( pj, nbpart, dt)
        push_x!( pj, nbpart, mesh, dt)

        f90_push_v!( pf, nbpart, dt)
        f90_push_x!( pf, nbpart, dimx, dimy, dt)

        @test all( pj .== pf)
    end

end
