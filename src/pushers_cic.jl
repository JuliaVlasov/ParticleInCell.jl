export push_x!, push_v!

function push_v!(p, kernel::CloudInCell, m::TwoDGrid, ex, ey, bz, dt)

    nbpart = size(p.array, 2)

    nx = m.nx
    ny = m.ny
    dx = m.dx
    dy = m.dy

    @threads for ipart = 1:nbpart

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

        ip1 = mod1(i+1, nx)
        jp1 = mod1(j+1, ny)

        e1 = a1 * ex[i, j] + a2 * ex[ip1, j] + a3 * ex[ip1, jp1] + a4 * ex[i, jp1]
        e2 = a1 * ey[i, j] + a2 * ey[ip1, j] + a3 * ey[ip1, jp1] + a4 * ey[i, jp1]
        b3 = a1 * bz[i, j] + a2 * bz[ip1, j] + a3 * bz[ip1, jp1] + a4 * bz[i, jp1]

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

function push_x!(p, mesh::TwoDGrid, dt::Float64)

    nbpart = size(p.array, 2)

    dimx, dimy = mesh.dimx, mesh.dimy

    @threads for i = 1:nbpart
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


function compute_rho(p, kernel::CloudInCell, m::TwoDGrid)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    rho = zeros(nx, ny)
    nbpart = size(p.array)[2]

    @inbounds for ipart = 1:nbpart
        xp = p.array[1, ipart] / dx
        yp = p.array[2, ipart] / dy

        i = trunc(Int, xp) + 1
        j = trunc(Int, yp) + 1

        dxp = xp - i + 1
        dyp = yp - j + 1

        a1 = (1 - dxp) * (1 - dyp)
        a2 = dxp * (1 - dyp)
        a3 = dxp * dyp
        a4 = (1 - dxp) * dyp

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        w = p.array[5, ipart]

        rho[i, j] += a1 * w
        rho[ip1, j] += a2 * w
        rho[ip1, jp1] += a3 * w
        rho[i, jp1] += a4 * w

    end

    rho ./= (dx * dy)

    return rho .- mean(rho)


end
