export compute_current!

function compute_current!(jx, jy, p, m)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    fill!(jx, 0)
    fill!(jy, 0)

    @inbounds for ipart = 1:p.nbpart

        xp = p.data[1, ipart] / dx
        yp = p.data[2, ipart] / dy

        i = trunc(Int, xp ) + 1
        j = trunc(Int, yp ) + 1

        dxp = i + 1 - xp
        dyp = j + 1 - yp
        dxq = 1 - dxp
        dyq = 1 - dyp

        a1 = dxp * dyp
        a2 = dxq * dyp
        a3 = dxq * dyq
        a4 = dxp * dyq

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        v1 = p.data[3, ipart]
        v2 = p.data[4, ipart]

        jx[i, j] += a1 * v1
        jy[i, j] += a1 * v2
        jx[ip1, j] += a2 * v1
        jy[ip1, j] += a2 * v2
        jx[ip1, jp1] += a3 * v1
        jy[ip1, jp1] += a3 * v2
        jx[i, jp1] += a4 * v1
        jy[i, jp1] += a4 * v2

    end

    jx .*= m.nx .* m.ny ./ p.nbpart ./ (dx * dy)
    jy .*= m.nx .* m.ny ./ p.nbpart ./ (dx * dy)

end
