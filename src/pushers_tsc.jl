export TriangularShapeCloud
  
f1_tsc( x ) = 0.75 - x^2
f2_tsc( x ) = 0.5 * ( 1.5 - x )^2

struct TriangularShapeCloud end

function push_v!(p, kernel::TriangularShapeCloud, m::TwoDGrid, ex, ey, bz, dt)

    nbpart = size(p.array, 2)

    nx = m.nx
    ny = m.ny
    dx = m.dx
    dy = m.dy

    @threads for ipart = 1:nbpart

        v1 = p.array[3, ipart]
        v2 = p.array[4, ipart]

        xp = p.array[1, ipart] / dx
        yp = p.array[2, ipart] / dy

        i = floor(Int, xp) + 1
        j = floor(Int, yp) + 1

        dpx = xp - i + 1
        dpy = yp - j + 1

        im1 = mod1(i-1,nx) 
        ip1 = mod1(i+1,nx) 
        jm1 = mod1(j-1,ny) 
        jp1 = mod1(j+1,ny) 

        cm1x = f2_tsc(1+dpx)
        cp1x = f2_tsc(1-dpx)
        cx   = f1_tsc(  dpx)
        cy   = f1_tsc(  dpy)
        cm1y = f2_tsc(1+dpy)
        cp1y = f2_tsc(1-dpy)

        e1 = 0.0
        e1 += cm1x * cm1y * ex[im1,jm1]   
        e1 += cm1x * cy   * ex[im1,j  ]   
        e1 += cm1x * cp1y * ex[im1,jp1]   
        e1 += cx   * cm1y * ex[i  ,jm1]   
        e1 += cx   * cy   * ex[i  ,j  ]   
        e1 += cx   * cp1y * ex[i  ,jp1]   
        e1 += cp1x * cm1y * ex[ip1,jm1]   
        e1 += cp1x * cy   * ex[ip1,j  ]   
        e1 += cp1x * cp1y * ex[ip1,jp1]   

        e2 = 0.0
        e2 += cm1x * cm1y * ey[im1,jm1]   
        e2 += cm1x * cy   * ey[im1,j  ]   
        e2 += cm1x * cp1y * ey[im1,jp1]   
        e2 += cx   * cm1y * ey[i  ,jm1]   
        e2 += cx   * cy   * ey[i  ,j  ]   
        e2 += cx   * cp1y * ey[i  ,jp1]   
        e2 += cp1x * cm1y * ey[ip1,jm1]   
        e2 += cp1x * cy   * ey[ip1,j  ]   
        e2 += cp1x * cp1y * ey[ip1,jp1]   

        b3 = 0.0
        b3 += cm1x * cm1y * bz[im1,jm1]   
        b3 += cm1x * cy   * bz[im1,j  ]   
        b3 += cm1x * cp1y * bz[im1,jp1]   
        b3 += cx   * cm1y * bz[i  ,jm1]   
        b3 += cx   * cy   * bz[i  ,j  ]   
        b3 += cx   * cp1y * bz[i  ,jp1]   
        b3 += cp1x * cm1y * bz[ip1,jm1]   
        b3 += cp1x * cy   * bz[ip1,j  ]   
        b3 += cp1x * cp1y * bz[ip1,jp1]   

        v1 += 0.5dt * e1
        v2 += 0.5dt * e2

        tantheta = 0.5dt * b3
        sintheta = 2 * tantheta / (1 + tantheta * tantheta)

        v1 += v2 * tantheta
        v2 += -v1 * sintheta
        v1 += v2 * tantheta

        p.array[3, ipart] = v1 + 0.5dt * e1
        p.array[4, ipart] = v2 + 0.5dt * e2

    end

end

function compute_rho( p , kernel::TriangularShapeCloud, mesh::TwoDGrid) 

    nbpart = size(p.array, 2)

    nx = mesh.nx
    ny = mesh.ny

    dx = mesh.dx
    dy = mesh.dy

    xmin, xmax = mesh.xmax, mesh.xmin
    ymin, ymax = mesh.ymax, mesh.ymin
    dimx = mesh.dimx
    dimy = mesh.dimy

    ρ = zeros(nx, ny)

    for ipart = 1:nbpart
    
        xp = p.array[1, ipart] / dx
        yp = p.array[2, ipart] / dy

        i = floor(Int, xp) + 1
        j = floor(Int, yp) + 1

        dpx = xp - i + 1
        dpy = yp - j + 1

        weight = p.array[5, ipart]
      
        im1 = mod1(i-1,nx)
        ip1 = mod1(i+1,nx)
        jm1 = mod1(j-1,ny)
        jp1 = mod1(j+1,ny)

        cm1x = f2_tsc(1+dpx)
        cp1x = f2_tsc(1-dpx)
        cx   = f1_tsc(  dpx)
        cy   = f1_tsc(  dpy)
        cm1y = f2_tsc(1+dpy)
        cp1y = f2_tsc(1-dpy)
      
        ρ[im1,jm1] += cm1x * cm1y * weight
        ρ[im1,j  ] += cm1x * cy   * weight
        ρ[im1,jp1] += cm1x * cp1y * weight
        ρ[i  ,jm1] += cx   * cm1y * weight
        ρ[i  ,j  ] += cx   * cy   * weight
        ρ[i  ,jp1] += cx   * cp1y * weight
        ρ[ip1,jm1] += cp1x * cm1y * weight
        ρ[ip1,j  ] += cp1x * cy   * weight
        ρ[ip1,jp1] += cp1x * cp1y * weight
      
    end
    
    ρ ./= (dx*dy)
    
    rho_total = sum(ρ) * dx * dy

    return ρ .- rho_total/dimx/dimy

end 

