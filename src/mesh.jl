struct Mesh

    nx
    ny
    x
    y
    hx
    hy
    hhx
    hhy

    function Mesh(dimx, nx, dimy, ny)

        x   = OffsetArray{Float64}(undef, -1:nx+1) 
        y   = OffsetArray{Float64}(undef, -1:ny+1)
        hx  = OffsetArray{Float64}(undef, -1:nx)
        hy  = OffsetArray{Float64}(undef, -1:ny)
        hhx = OffsetArray{Float64}(undef,  0:nx)
        hhy = OffsetArray{Float64}(undef,  0:ny)
        
        dx = dimx / nx
        dy = dimy / ny
        
        x[0] = 0
        y[0] = 0
        
        for i=1:nx
            x[i] = (i*dx) *(i*dx+1)/(1+dimx)
        end
        
        for j=1:ny
            y[j] = (j*dy) *(j*dy+1)/(1+dimy)
        end
        
        for i=0:nx-1
            hx[i] = x[i+1]-x[i]
        end
        
        for j=0:ny-1
            hy[j] = y[j+1]-y[j]
        end
        
        hx[nx] = hx[0]  
        hx[-1] = hx[nx-1]
        hy[ny] = hy[0]
        hy[-1] = hy[ny-1]
        
        x[-1]   = x[0] - hx[nx-1]  # points utiles pour le cas period
        x[nx+1] = x[nx] + hx[0]
        y[-1]   = y[0] - hy[ny-1]
        y[ny+1] = y[ny] + hy[0]
        
        hhx[0]  =  0.5 * ( hx[0] + hx[nx-1] ) 
        hhx[nx] =  0.5 * ( hx[0] + hx[nx-1] ) 
        for i=1:nx-1
           hhx[i] = 0.5 * ( hx[i] + hx[i-1] )
        end
        
        hhy[0]  = 0.5 * ( hy[0] + hy[ny-1] ) 
        hhy[ny] = 0.5 * ( hy[0] + hy[ny-1] ) 
        
        for j=1:ny-1
           hhy[j] = 0.5 * ( hy[j] + hy[j-1] )
        end

        new( nx, ny, x, y, hx, hy, hhx, hhy)

    end

end
