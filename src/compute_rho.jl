using Statistics

export compute_rho

function compute_rho(p, m)

    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    rho = zeros(nx, ny)
    nbpart = size(p)[2]

    for ipart = 1:nbpart
        xp = p[1, ipart] / dx
        yp = p[2, ipart] / dy

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

        rho[i, j] += a1
        rho[ip1, j] += a2
        rho[ip1, jp1] += a3
        rho[i, jp1] += a4

    end

    rho = rho .* (nx * ny) ./ nbpart

    return rho .- mean(rho)

end
