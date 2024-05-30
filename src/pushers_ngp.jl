export NearestGridPoint
  
struct NearestGridPoint end

function push_v!(p, kernel::NearestGridPoint, m::TwoDGrid, ex, ey, bz, dt)

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

        i = mod1(round(Int, xp) + 1, nx)
        j = mod1(round(Int, yp) + 1, ny)

        e1 = ex[i,j]   
        e2 = ey[i,j]   
        b3 = bz[i,j]   

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

function compute_rho( p , kernel::NearestGridPoint, mesh::TwoDGrid) 

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

        i = mod1(round(Int, xp) + 1, nx)
        j = mod1(round(Int, yp) + 1, ny)

        weight = p.array[5, ipart]
      
        ρ[i,j] += weight
      
    end
    
    ρ ./= (dx*dy)
    
    rho_total = sum(ρ) * dx * dy

    return ρ .- rho_total/dimx/dimy

end 
