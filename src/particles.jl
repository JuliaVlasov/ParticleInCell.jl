import Base.Threads

export push_v!

function push_v!(p, kernel:: CloudInCell, m :: TwoDGrid, dt)

    nbpart = size(p.array)[2]

    dx = m.dx
    dy = m.dy

    Threads.@threads for ipart = 1:nbpart

        v1 = p.array[3, ipart]
        v2 = p.array[4, ipart]

        xp = p.array[1, ipart] / dx
        yp = p.array[2, ipart] / dy

        i = floor(Int, xp) + 1
        j = floor(Int, yp) + 1

        dxp = xp - i + 1
        dyp = yp - j + 1

        a1 = (1 - dxp) * (1 - dyp)
        a2 = dxp * (1 - dyp)
        a3 = dxp * dyp
        a4 = (1 - dxp) * dyp

        e1 = a1 * m.ex[i, j] + a2 * m.ex[i+1, j] + a3 * m.ex[i+1, j+1] + a4 * m.ex[i, j+1]
        e2 = a1 * m.ey[i, j] + a2 * m.ey[i+1, j] + a3 * m.ey[i+1, j+1] + a4 * m.ey[i, j+1]
        b3 = a1 * m.bz[i, j] + a2 * m.bz[i+1, j] + a3 * m.bz[i+1, j+1] + a4 * m.bz[i, j+1]

        v1 += 0.5dt * e1
        v2 += 0.5dt * e2

        tantheta = 0.5dt * b3
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 += v2 * tantheta
        v2 += -v1 * sintheta
        v1 += v2 * tantheta

        p.array[3, ipart] = v1 + 0.5dt * e1
        p.array[4, ipart] = v2 + 0.5dt * e2

    end

end

export push_x!

function push_x!(p, mesh::TwoDGrid, dt::Float64)

    nbpart = size(p.array)[2]

    dimx, dimy = mesh.dimx, mesh.dimy

    Threads.@threads for i = 1:nbpart
        p1 = p.array[1, i] + dt * p.array[3, i]
        p2 = p.array[2, i] + dt * p.array[4, i]
        p1 > dimx && (p1 -= dimx)
        p2 > dimy && (p2 -= dimy)
        p1 < 0.0 && (p1 += dimx)
        p2 < 0.0 && (p2 += dimy)
        p.array[1, i] = p1
        p.array[2, i] = p2
    end

end
