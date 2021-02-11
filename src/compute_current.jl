export compute_current!

function compute_current!(jx, jy, p, m)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    fill!(jx, 0)
    fill!(jy, 0)

    @inbounds for ipart = 1:p.nbpart
        i = p.cell[1,ipart]
        j = p.cell[2,ipart]

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        xp = p.pos[1, ipart]
        yp = p.pos[2, ipart]

        a1 = (m.x[i+1] - xp) * (m.y[j+1] - yp)
        a2 = (xp - m.x[i]) * (m.y[j+1] - yp)
        a3 = (xp - m.x[i]) * (yp - m.y[j])
        a4 = (m.x[i+1] - xp) * (yp - m.y[j])

        v1 = p.vit[1, ipart]

        jx[i, j] += a1 * v1
        jx[ip1, j] += a2 * v1
        jx[ip1, jp1] += a3 * v1
        jx[i, jp1] += a4 * v1

        v2 = p.vit[2, ipart]

        jy[i, j] += a1 * v2
        jy[ip1, j] += a2 * v2
        jy[ip1, jp1] += a3 * v2
        jy[i, jp1] += a4 * v2

    end

    jx .= jx .* m.nx .* m.ny ./ p.nbpart ./ (dx * dy)
    jy .= jy .* m.nx .* m.ny ./ p.nbpart ./ (dx * dy)

end
