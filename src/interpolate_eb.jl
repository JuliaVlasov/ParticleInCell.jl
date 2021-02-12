export interpol_eb!

function interpol_eb!(eb, p::Particles, m::Mesh)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    @inbounds for ip = 1:p.nbpart

        xp = p.pos[1, ip]
        yp = p.pos[2, ip]

        i = trunc(Int, xp / dx) + 1
        j = trunc(Int, yp / dy) + 1

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        dxp = i*dx - xp
        dyp = j*dy - yp
        dxq = dx - dxp
        dyq = dy - dyp

        a1 = dxp * dyp
        a2 = dxq * dyp
        a3 = dxq * dyq
        a4 = dxp * dyq

        p.ebp[1,ip] = a1 * eb[1,i,j] + a2 * eb[1,ip1,j] + a3 * eb[1,ip1,jp1] + a4 * eb[1,i,jp1]
        p.ebp[2,ip] = a1 * eb[2,i,j] + a2 * eb[2,ip1,j] + a3 * eb[2,ip1,jp1] + a4 * eb[2,i,jp1]
        p.ebp[3,ip] = a1 * eb[3,i,j] + a2 * eb[3,ip1,j] + a3 * eb[3,ip1,jp1] + a4 * eb[3,i,jp1]

    end

    p.ebp ./= (m.dx * m.dy)

end
