export interpolate_cic!

function interpolate_cic!(epx, epy, mesh, p, ex, ey)

    dx = mesh.dx
    dy = mesh.dy

    for m = eachindex(epx, epy)
    
       xp = p.array[1, m] / dx
       yp = p.array[2, m] / dy

       i = trunc(Int, xp) + 1
       j = trunc(Int, yp) + 1

       dxp = xp - i + 1
       dyp = yp - j + 1

       a1 = (1 - dxp) * (1 - dyp)
       a2 = dxp * (1 - dyp)
       a3 = dxp * dyp
       a4 = (1 - dxp) * dyp

       ip1 = mod1(i + 1, mesh.nx)
       jp1 = mod1(j + 1, mesh.ny)

       epx[m] = a1 * ex[i,  j] + a2 * ex[ip1, j] + a3 * ex[ip1,  jp1] + a4 * ex[i, jp1]
       epy[m] = a1 * ey[i,  j] + a2 * ey[ip1, j] + a3 * ey[ip1,  jp1] + a4 * ey[i, jp1]
                
    
    end

end
