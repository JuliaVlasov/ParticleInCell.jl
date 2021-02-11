export interpol_eb!

function interpol_eb!(ex, ey, bz, p::Particles, m::Mesh)

    dum = 1 / (m.dx * m.dy)

    @inbounds for ipart = 1:p.nbpart

        i = p.cell[1,ipart]
        j = p.cell[2,ipart]

        xp = p.pos[1, ipart]
        yp = p.pos[2, ipart]

        ip1 = mod1(i + 1, m.nx)
        jp1 = mod1(j + 1, m.ny)

        a1 = (m.x[i+1] - xp) * (m.y[j+1] - yp) * dum
        a2 = (xp - m.x[i]) * (m.y[j+1] - yp) * dum
        a3 = (xp - m.x[i]) * (yp - m.y[j]) * dum
        a4 = (m.x[i+1] - xp) * (yp - m.y[j]) * dum

        p.epx[ipart] = a1 * ex[i, j] + a2 * ex[ip1, j] + a3 * ex[ip1, jp1] + a4 * ex[i, jp1]
        p.epy[ipart] = a1 * ey[i, j] + a2 * ey[ip1, j] + a3 * ey[ip1, jp1] + a4 * ey[i, jp1]
        p.bpz[ipart] = a1 * bz[i, j] + a2 * bz[ip1, j] + a3 * bz[ip1, jp1] + a4 * bz[i, jp1]

    end

end
