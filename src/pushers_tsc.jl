export TriangularShapeCloud
  
function f_tsc( x ) 

     if x <= 0.5 
         0.75 - x * x
     elseif 0.5 < x <= 1.5 
         0.5 * (1.5 - x)^2
     else
         0.0
     end

end

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

        im2 = mod1(i-2,nx) 
        im1 = mod1(i-1,nx) 
        ip1 = mod1(i+1,nx) 
        ip2 = mod1(i+2,nx) 

        jm2 = mod1(j-2,ny) 
        jm1 = mod1(j-1,ny) 
        jp1 = mod1(j+1,ny) 
        jp2 = mod1(j+2,ny) 

        cm2x = f_m4(2+dpx)
        cp2x = f_m4(2-dpx)
        cm1x = f_m4(1+dpx)
        cp1x = f_m4(1-dpx)
        cx   = f_m4(  dpx)
        cy   = f_m4(  dpy)
        cm2y = f_m4(2+dpy)
        cp2y = f_m4(2-dpy)
        cm1y = f_m4(1+dpy)
        cp1y = f_m4(1-dpy)

        e1 = 0.0
        e1 += cm2x * cm2y * ex[im2,jm2]   
        e1 += cm2x * cm1y * ex[im2,jm1]   
        e1 += cm2x * cy   * ex[im2,j  ]   
        e1 += cm2x * cp1y * ex[im2,jp1]   
        e1 += cm2x * cp2y * ex[im2,jp2]   
        e1 += cm1x * cm2y * ex[im1,jm2]   
        e1 += cm1x * cm1y * ex[im1,jm1]   
        e1 += cm1x * cy   * ex[im1,j  ]   
        e1 += cm1x * cp1y * ex[im1,jp1]   
        e1 += cm1x * cp2y * ex[im1,jp2]   
        e1 += cx   * cm2y * ex[i  ,jm2]   
        e1 += cx   * cm1y * ex[i  ,jm1]   
        e1 += cx   * cy   * ex[i  ,j  ]   
        e1 += cx   * cp1y * ex[i  ,jp1]   
        e1 += cx   * cp2y * ex[i  ,jp2]   
        e1 += cp1x * cm2y * ex[ip1,jm2]   
        e1 += cp1x * cm1y * ex[ip1,jm1]   
        e1 += cp1x * cy   * ex[ip1,j  ]   
        e1 += cp1x * cp1y * ex[ip1,jp1]   
        e1 += cp1x * cp2y * ex[ip1,jp2]   
        e1 += cp2x * cm2y * ex[ip2,jm2]   
        e1 += cp2x * cm1y * ex[ip2,jm1]   
        e1 += cp2x * cy   * ex[ip2,j  ]   
        e1 += cp2x * cp1y * ex[ip2,jp1]   
        e1 += cp2x * cp2y * ex[ip2,jp2]   

        e2 = 0.0
        e2 += cm2x * cm2y * ey[im2,jm2]   
        e2 += cm2x * cm1y * ey[im2,jm1]   
        e2 += cm2x * cy   * ey[im2,j  ]   
        e2 += cm2x * cp1y * ey[im2,jp1]   
        e2 += cm2x * cp2y * ey[im2,jp2]   
        e2 += cm1x * cm2y * ey[im1,jm2]   
        e2 += cm1x * cm1y * ey[im1,jm1]   
        e2 += cm1x * cy   * ey[im1,j  ]   
        e2 += cm1x * cp1y * ey[im1,jp1]   
        e2 += cm1x * cp2y * ey[im1,jp2]   
        e2 += cx   * cm2y * ey[i  ,jm2]   
        e2 += cx   * cm1y * ey[i  ,jm1]   
        e2 += cx   * cy   * ey[i  ,j  ]   
        e2 += cx   * cp1y * ey[i  ,jp1]   
        e2 += cx   * cp2y * ey[i  ,jp2]   
        e2 += cp1x * cm2y * ey[ip1,jm2]   
        e2 += cp1x * cm1y * ey[ip1,jm1]   
        e2 += cp1x * cy   * ey[ip1,j  ]   
        e2 += cp1x * cp1y * ey[ip1,jp1]   
        e2 += cp1x * cp2y * ey[ip1,jp2]   
        e2 += cp2x * cm2y * ey[ip2,jm2]   
        e2 += cp2x * cm1y * ey[ip2,jm1]   
        e2 += cp2x * cy   * ey[ip2,j  ]   
        e2 += cp2x * cp1y * ey[ip2,jp1]   
        e2 += cp2x * cp2y * ey[ip2,jp2]   

        b3 = 0.0
        b3 += cm2x * cm2y * bz[im2,jm2]   
        b3 += cm2x * cm1y * bz[im2,jm1]   
        b3 += cm2x * cy   * bz[im2,j  ]   
        b3 += cm2x * cp1y * bz[im2,jp1]   
        b3 += cm2x * cp2y * bz[im2,jp2]   
        b3 += cm1x * cm2y * bz[im1,jm2]   
        b3 += cm1x * cm1y * bz[im1,jm1]   
        b3 += cm1x * cy   * bz[im1,j  ]   
        b3 += cm1x * cp1y * bz[im1,jp1]   
        b3 += cm1x * cp2y * bz[im1,jp2]   
        b3 += cx   * cm2y * bz[i  ,jm2]   
        b3 += cx   * cm1y * bz[i  ,jm1]   
        b3 += cx   * cy   * bz[i  ,j  ]   
        b3 += cx   * cp1y * bz[i  ,jp1]   
        b3 += cx   * cp2y * bz[i  ,jp2]   
        b3 += cp1x * cm2y * bz[ip1,jm2]   
        b3 += cp1x * cm1y * bz[ip1,jm1]   
        b3 += cp1x * cy   * bz[ip1,j  ]   
        b3 += cp1x * cp1y * bz[ip1,jp1]   
        b3 += cp1x * cp2y * bz[ip1,jp2]   
        b3 += cp2x * cm2y * bz[ip2,jm2]   
        b3 += cp2x * cm1y * bz[ip2,jm1]   
        b3 += cp2x * cy   * bz[ip2,j  ]   
        b3 += cp2x * cp1y * bz[ip2,jp1]   
        b3 += cp2x * cp2y * bz[ip2,jp2]   

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
      
        im2 = mod1(i-2,nx)
        im1 = mod1(i-1,nx)
        ip1 = mod1(i+1,nx)
        ip2 = mod1(i+2,nx)
        jm2 = mod1(j-2,ny)
        jm1 = mod1(j-1,ny)
        jp1 = mod1(j+1,ny)
        jp2 = mod1(j+2,ny)

        cm2x = f_tsc(2+dpx)
        cp2x = f_tsc(2-dpx)
        cm1x = f_tsc(1+dpx)
        cp1x = f_tsc(1-dpx)
        cx   = f_tsc(dpx)
        cy   = f_tsc(dpy)
        cp1y = f_tsc(1-dpy)
        cm1y = f_tsc(1+dpy)
        cp2y = f_tsc(2-dpy)
        cm2y = f_tsc(2+dpy)
      
        ρ[im2,jm2] += cm2x * cm2y * weight
        ρ[im2,jm1] += cm2x * cm1y * weight
        ρ[im2,j  ] += cm2x * cy   * weight
        ρ[im2,jp1] += cm2x * cp1y * weight
        ρ[im2,jp2] += cm2x * cp2y * weight
      
        ρ[im1,jm2] += cm1x * cm2y * weight
        ρ[im1,jm1] += cm1x * cm1y * weight
        ρ[im1,j  ] += cm1x * cy   * weight
        ρ[im1,jp1] += cm1x * cp1y * weight
        ρ[im1,jp2] += cm1x * cp2y * weight

        ρ[i  ,jm2] += cx   * cm2y * weight
        ρ[i  ,jm1] += cx   * cm1y * weight
        ρ[i  ,j  ] += cx   * cy   * weight
        ρ[i  ,jp1] += cx   * cp1y * weight
        ρ[i  ,jp2] += cx   * cp2y * weight
      
        ρ[ip1,jm2] += cp1x * cm2y * weight
        ρ[ip1,jm1] += cp1x * cm1y * weight
        ρ[ip1,j  ] += cp1x * cy   * weight
        ρ[ip1,jp1] += cp1x * cp1y * weight
        ρ[ip1,jp2] += cp1x * cp2y * weight
      
        ρ[ip2,jm2] += cp2x * cm2y * weight
        ρ[ip2,jm1] += cp2x * cm1y * weight
        ρ[ip2,j  ] += cp2x * cy   * weight
        ρ[ip2,jp1] += cp2x * cp1y * weight
        ρ[ip2,jp2] += cp2x * cp2y * weight
      
    end
    
    ρ ./= (dx*dy)
    
    rho_total = sum(ρ) * dx * dy

    return ρ .- rho_total/dimx/dimy

end 
