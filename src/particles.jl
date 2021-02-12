export Particles

struct Particles

    nbpart::Int
    data::Array{Float64,2}

    function Particles(nbpart)

        data = zeros(7, nbpart)
        new(nbpart, data)

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
        pg.data[1, i] = newton(r2)
        pg.data[2, i] = r3
        pg.data[3, i] = v * cos(θ)
        pg.data[4, i] = v * sin(θ)
    end

end

export update_cells!

function update_cells!(p, m)

    p.data[1,:] .= mod.(p.data[1,:], m.dimx)
    p.data[2,:] .= mod.(p.data[2,:], m.dimy)

end


export push_v!

function push_v!(p, dt)

    for ipart = 1:p.nbpart

        v1 = p.data[3, ipart]
        v2 = p.data[4, ipart]
        e1 = p.data[5, ipart]
        e2 = p.data[6, ipart]
        b3 = p.data[7, ipart]

        v1 += 0.5dt * e1
        v2 += 0.5dt * e2

        tantheta = 0.5dt * b3
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 += v2 * tantheta
        v2 += - v1 * sintheta
        v1 += v2 * tantheta

        p.data[3, ipart] = v1 + 0.5dt * e1
        p.data[4, ipart] = v2 + 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh, dt)

    @simd for ipart in 1:p.nbpart
        p.data[1,ipart] += p.data[3,ipart] * dt
        p.data[2,ipart] += p.data[4,ipart] * dt
    end

    update_cells!(p, mesh)

end
