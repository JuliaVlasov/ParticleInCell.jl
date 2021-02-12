export FDTD

struct FDTD

    ex
    ey
    bz
    jx
    jy

    function FDTD( m )

        nx, ny = m.nx, m.ny
        ex = zeros(nx,ny)
        ey = zeros(nx,ny)
        bz = zeros(nx,ny)
        jx = zeros(nx,ny)
        jy = zeros(nx,ny)

        new( ex, ey, bz, jx, jy )

    end

end

function faraday!( eb, fdtd :: FDTD, m, dt )

   nx, ny = m.nx, m.ny
   dx, dy = m.dx, m.dy

   for j=1:ny, i=1:nx
      dex_dy  = (fdtd.ex[i,mod1(j+1,ny)]-fdtd.ex[i,j]) / dy
      dey_dx  = (fdtd.ey[mod1(i+1,nx),j]-fdtd.ey[i,j]) / dx
      fdtd.bz[i,j] += dt * (dex_dy - dey_dx)
   end

   for j=1:ny, i=1:nx
      eb[3,i,j] = ( fdtd.bz[mod1(i-1,nx), mod1(j-1,ny)] + fdtd.bz[i, mod1(j-1,ny)] 
                  + fdtd.bz[mod1(i-1,nx), j] + fdtd.bz[i, j]) / 4
   end
   
end 

function ampere_maxwell!( ex, ey, fdtd :: FDTD, m, dt)

   nx, ny = m.nx, m.ny
   dx, dy = m.dx, m.dy

   for i=1:nx, j=1:ny
      dbz_dy = (fdtd.bz[i,j]-fdtd.bz[i,mod1(j-1,ny)]) / dy
      fdtd.ex[i,j] += dt * dbz_dy - dt * fdtd.jx[i,j]
   end

   for i=1:nx, j=1:ny
      dbz_dx = (fdtd.bz[i,j]-fdtd.bz[mod1(i-1,nx),j]) / dx
      fdtd.ey[i,j] -= dt * dbz_dx - dt * fdtd.jy[i,j]
   end

   for i=1:nx, j=1:ny
      ex[i,j] = 0.5*( fdtd.ex[mod1(i-1,nx),j] + fdtd.ex[i,j] )
      ey[i,j] = 0.5*( fdtd.ey[i,mod1(j-1,ny)] + fdtd.ey[i,j] )
   end

end 



function compute_current!( jx, jy, fdtd :: FDTD, p, m )


    nx, ny = m.nx, m.ny
    dx, dy = m.dx, m.dy
    
    fill!(jx, 0)
    fill!(jy, 0)
    
    for ipart=1:p.nbpart
    
       xp = p.pos[1,ipart]
       yp = p.pos[2,ipart]

       i = trunc(Int, xp / dx) + 1
       j = trunc(Int, yp / dy) + 1
    
       ip1 = mod1(i+1, nx)
       jp1 = mod1(j+1, ny)
    
       dxp = i*dx - xp
       dyp = j*dy - yp
       dxq = dx - dxp
       dyq = dy - dyp

       a1 = dxp * dyp
       a2 = dxq * dyp
       a3 = dxq * dyq
       a4 = dxp * dyq
    
       w = p.vit[1,ipart]
    
       jx[i,j]     +=  a1*w  
       jx[ip1,j]   +=  a2*w 
       jx[ip1,jp1] +=  a3*w 
       jx[i,jp1]   +=  a4*w 
    
       w = p.vit[2,ipart]
    
       jy[i,j]     += a1*w  
       jy[ip1,j]   += a2*w 
       jy[ip1,jp1] += a3*w 
       jy[i,jp1]   += a4*w 
    
    end

    jx .*= (nx * ny) / p.nbpart / (dx * dy)
    jy .*= (nx * ny) / p.nbpart / (dx * dy)

    
    for i=1:nx, j=1:ny
       fdtd.jx[i,j] = 0.5 * (jx[i,j]+jx[mod1(i+1,nx),j])
    end
    
    for i=1:nx, j=1:ny
       fdtd.jy[i,j] = 0.5 * (jy[i,j]+jy[i,mod1(j+1,ny)])
    end

end 

