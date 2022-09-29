using Sobol

export landau_sampling

function landau_sampling( mesh :: Mesh, nbpart :: Int64 )

    kx     = 0.5
    alpha  = 0.05

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx = mesh.xmax - mesh.xmin
    dimy = mesh.ymax - mesh.ymin

    weight = (dimx * dimy) / nbpart

    particles = Particles(nbpart, weight )

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1-x0) > 1e-12)
            p = x0 + alpha * sin( kx * x0) / kx 
            f = 1 + alpha * cos( kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end
    
    s = Sobol.SobolSeq(3)
    nbpart = pg.n_particles

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2, r3  = Sobol.next!(s)
        θ = r1 * 2π
        particles.x[1,i] = xmin + newton(r2) * ( xmax - xmin )
        particles.x[2,i] = ymin + r3 * ( ymax - ymin )
        particles.v[1,i] = v * cos(θ)
        particles.v[2,i] = v * sin(θ)

    end

    particles

end
