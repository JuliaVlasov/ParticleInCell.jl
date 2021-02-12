using Statistics

export compute_rho

function compute_rho(p, m)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    rho = zeros(nx, ny)

    for ipart = 1:p.nbpart

        xp = p.data[1, ipart] / dx
        yp = p.data[2, ipart] / dy

        i = trunc(Int, xp ) + 1
        j = trunc(Int, yp ) + 1

        dxp = i - xp
        dyp = j - yp

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        a1 = dxp * dyp
        a2 = (1 - dxp) * dyp
        a3 = (1 - dxp) * (1 - dyp)
        a4 = dxp * (1 - dyp)

        rho[i, j] += a1
        rho[ip1, j] += a2
        rho[ip1, jp1] += a3
        rho[i, jp1] += a4

    end

    rho = rho .* (nx * ny) ./ p.nbpart 

    return rho .- mean(rho)

end
