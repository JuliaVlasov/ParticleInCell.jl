export landau_sampling!

function landau_sampling!(particles, alpha, kx)

    nbpart = size(particles)[2]

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1 - x0) > 1e-12)
            p = x0 + alpha * sin(kx * x0) / kx
            f = 1 + alpha * cos(kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end

    s = Sobol.SobolSeq(3)

    for i = 1:nbpart
        v = sqrt(-2 * log((i - 0.5) / nbpart))
        r1, r2, r3 = Sobol.next!(s)
        θ = r1 * 2π
        particles[1, i] = newton(r2)
        particles[2, i] = r3
        particles[3, i] = v * cos(θ)
        particles[4, i] = v * sin(θ)
    end

end

export update_cells!

export push_v!

function push_v!(p, dt)

    nbpart = size(p)[2]

    for ipart = 1:nbpart

        v1 = p[3, ipart]
        v2 = p[4, ipart]
        e1 = p[5, ipart]
        e2 = p[6, ipart]
        b3 = p[7, ipart]

        v1 += 0.5dt * e1
        v2 += 0.5dt * e2

        tantheta = 0.5dt * b3
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 += v2 * tantheta
        v2 += - v1 * sintheta
        v1 += v2 * tantheta

        p[3, ipart] = v1 + 0.5dt * e1
        p[4, ipart] = v2 + 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh :: Mesh, dt :: Float64)

    nbpart = size(p)[2]

    dimx, dimy = mesh.dimx, mesh.dimy

    p[1,:] .= mod.(view(p,1,:) .+ dt .* view(p,3,:), dimx)
    p[2,:] .= mod.(view(p,2,:) .+ dt .* view(p,4,:), dimy)

end
