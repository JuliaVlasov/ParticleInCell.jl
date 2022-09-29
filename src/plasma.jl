export plasma

function plasma( mesh :: Mesh, nbpart :: Int64 )

    kx     = 0.5
    alpha  = 0.05

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx = mesh.xmax - mesh.xmin
    dimy = mesh.ymax - mesh.ymin

    weight = (dimx * dimy) / nbpart

    particles = Particles(nbpart, weight )

    k = 1
    while (k<=nbpart)

        xi   = rand() * dimx
        yi   = rand() * dimy
        zi   = (2.0 + alpha) * rand()
        temm = 1.0 + sin(yi) + alpha * cos(kx*xi)

        if (temm >= zi)
            particles.x[1,k] = xi
            particles.x[2,k] = yi
            k = k + 1
        end

    end

    k = 1
    while (k<=nbpart)

        xi   = (rand()-0.5)*10
        yi   = (rand()-0.5)*10
        zi   = rand()
        temm = ( exp(-((xi-2)^2 + yi^2)/2)
               + exp(-((xi+2)^2 + yi^2)/2))/2

        if (temm >= zi)
            particles.v[1,k] = xi
            particles.v[2,k] = yi
            k = k + 1
        end

    end

    particles

end
