struct Particle

   cell :: Int32
   dx   :: Float32
   v    :: Float64
   w    :: Float32

end

export Particles

struct Particles

    nbpart
    pos
    cell
    vit
    epx
    epy
    bpz
    p

    function Particles( nbpart )

        pos = zeros(nbpart,2)
        cell = zeros(Int32, nbpart,2)
        vit = zeros(nbpart,2)
        epx = zeros(nbpart)
        epy = zeros(nbpart)
        bpz = zeros(nbpart)
        p = zeros(nbpart)

        new( nbpart, pos, cell, vit, epx, epy, bpz, p )

    end 

end

export landau_sampling!

function landau_sampling!( pg, alpha, kx )

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
    nbpart = pg.nbpart

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2, r3 = Sobol.next!(s)
        θ = r1 * 2π
        pg.pos[i,1] = newton(r2)
        pg.pos[i,2] = r3 
        pg.vit[i,1] = v * cos(θ)
        pg.vit[i,2] = v * sin(θ)
        pg.p[i] = 1 / nbpart
    end

end

export update_cells!

function update_cells!( p, m )

    for i in 1:p.nbpart

       p.pos[i,1] = mod(p.pos[i,1], m.dimx)
       p.pos[i,2] = mod(p.pos[i,2], m.dimy)

       p.cell[i,1] = trunc(Int, p.pos[i,1] / m.dx) + 1
       p.cell[i,2] = trunc(Int, p.pos[i,2] / m.dy) + 1

    end

end
