export interpol_eb!

function interpol_eb!(ex, ey, bz, p::Particles, m::Mesh)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy

    @inbounds for ip = 1:p.nbpart

        xp = p.data[1, ip] / dx
        yp = p.data[2, ip] / dy

        i = trunc(Int, xp ) + 1
        j = trunc(Int, yp ) + 1

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        dxp = i - 1 - xp
        dyp = j - 1 - yp
        dxq = 1 - dxp
        dyq = 1 - dyp

        a1 = dxp * dyp
        a2 = dxq * dyp
        a3 = dxq * dyq
        a4 = dxp * dyq

        p.ebp[1,ip] = a1 * ex[i,j] + a2 * ex[ip1,j] + a3 * ex[ip1,jp1] + a4 * ex[i,jp1]
        p.ebp[2,ip] = a1 * ey[i,j] + a2 * ey[ip1,j] + a3 * ey[ip1,jp1] + a4 * ey[i,jp1]
        p.ebp[3,ip] = a1 * bz[i,j] + a2 * bz[ip1,j] + a3 * bz[ip1,jp1] + a4 * bz[i,jp1]

    end


end
