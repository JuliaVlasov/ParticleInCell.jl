export compute_current!

function compute_current!(jxy, p, m)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    fill!(jxy, 0)

    @inbounds for ipart = 1:p.nbpart


        xp = p.pos[1, ipart]
        yp = p.pos[2, ipart]

        i = trunc(Int, xp / m.dx) + 1
        j = trunc(Int, yp / m.dy) + 1

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        a1 = (m.x[i+1] - xp) * (m.y[j+1] - yp)
        a2 = (xp - m.x[i]) * (m.y[j+1] - yp)
        a3 = (xp - m.x[i]) * (yp - m.y[j])
        a4 = (m.x[i+1] - xp) * (yp - m.y[j])

        v1 = p.vit[1, ipart]

        jxy[1, i, j] += a1 * v1
        jxy[1, ip1, j] += a2 * v1
        jxy[1, ip1, jp1] += a3 * v1
        jxy[1, i, jp1] += a4 * v1

        v2 = p.vit[2, ipart]

        jxy[2, i, j] += a1 * v2
        jxy[2, ip1, j] += a2 * v2
        jxy[2, ip1, jp1] += a3 * v2
        jxy[2, i, jp1] += a4 * v2

    end

    jx .*= m.nx .* m.ny ./ p.nbpart ./ (dx * dy)
    jy .*= m.nx .* m.ny ./ p.nbpart ./ (dx * dy)

end
