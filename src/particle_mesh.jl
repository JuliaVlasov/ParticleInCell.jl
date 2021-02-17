export interpolation!

function interpolation!(p :: Array{Float64,2}, m :: Mesh)

    nbpart = size(p)[2]
    dx = m.dx
    dy = m.dy

    @inbounds for ipart = 1:nbpart

        xp = p[1, ipart] / dx
        yp = p[2, ipart] / dy

        i = floor(Int, xp ) + 1
        j = floor(Int, yp ) + 1

        dxp = xp - i + 1
        dyp = yp - j + 1
    
        a1 = (1-dxp) * (1-dyp)
        a2 = dxp * (1-dyp)
        a3 = dxp * dyp
        a4 = (1-dxp) * dyp

        p[5,ipart] = a1 * m.ex[ i, j] + a2 * m.ex[ i+1, j] + a3 * m.ex[ i+1, j+1] + a4 * m.ex[ i, j+1]
        p[6,ipart] = a1 * m.ey[ i, j] + a2 * m.ey[ i+1, j] + a3 * m.ey[ i+1, j+1] + a4 * m.ey[ i, j+1]
        p[7,ipart] = a1 * m.bz[ i, j] + a2 * m.bz[ i+1, j] + a3 * m.bz[ i+1, j+1] + a4 * m.bz[ i, j+1]

    end


end

export compute_current!

function compute_current!( m :: Mesh, p )

    nbpart = size(p)[2]
    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    
    fill!(m.jx, 0)
    fill!(m.jy, 0)

    factor = ( nx * ny ) / nbpart
    
    for ipart=1:nbpart
    
       xp = p[1,ipart] / dx
       yp = p[2,ipart] / dy

       i = floor(Int, xp) + 1
       j = floor(Int, yp) + 1
    
       dxp = xp - i + 1
       dyp = yp - j + 1
    
       a1 = (1-dxp) * (1-dyp)
       a2 = dxp * (1-dyp)
       a3 = dxp * dyp
       a4 = (1-dxp) * dyp

       w1 = p[3,ipart] * factor
       w2 = p[4,ipart] * factor
    
       m.jx[i,j]     += a1*w1
       m.jy[i,j]     += a1*w2  

       m.jx[i+1,j]   += a2*w1 
       m.jy[i+1,j]   += a2*w2 

       m.jx[i+1,j+1] += a3*w1 
       m.jy[i+1,j+1] += a3*w2 

       m.jx[i,j+1]   += a4*w1 
       m.jy[i,j+1]   += a4*w2 
    
    end

    for i=1:nx+1
      m.jx[i,1] += m.jx[i,ny+1]
      m.jx[i,ny+1] = m.jx[i,1]
    end
    for j=1:ny+1
      m.jx[1,j] += m.jx[nx+1,j]
      m.jx[nx+1,j] = m.jx[1,j]
    end

end 
