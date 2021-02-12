export Particles

struct Particles

    nbpart::Any
    pos::Any
    vit::Any
    ebp::Any

    function Particles(nbpart)

        pos = zeros(2, nbpart)
        vit = zeros(2, nbpart)
        ebp = zeros(3, nbpart)

        new(nbpart, pos, vit, ebp)

    end

end

export landau_sampling!

function landau_sampling!(pg, alpha, kx)

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
    nbpart = pg.nbpart

    for i = 1:nbpart
        v = sqrt(-2 * log((i - 0.5) / nbpart))
        r1, r2, r3 = Sobol.next!(s)
        θ = r1 * 2π
        pg.pos[1, i] = newton(r2)
        pg.pos[2, i] = r3
        pg.vit[1, i] = v * cos(θ)
        pg.vit[2, i] = v * sin(θ)
    end

end

export update_cells!

function update_cells!(p, m)

    p.pos[1,:] .= mod.(p.pos[1,:], m.dimx)
    p.pos[2,:] .= mod.(p.pos[2,:], m.dimy)

end
