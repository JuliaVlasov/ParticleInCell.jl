export FDTD

struct FDTD

    ex
    ey
    bz
    jx
    jy

    function FDTD( mesh )

        nx, ny = mesh.nx, mesh.ny
        ex = zeros(nx,ny)
        ey = zeros(nx,ny)
        bz = zeros(nx,ny)
        jx = zeros(nx,ny)
        jy = zeros(nx,ny)

        new( ex, ey, bz, jx, jy )

    end

end

function faraday!( bz, fdtd :: FDTD, mesh, dt )

   nx, ny = mesh.nx, mesh.ny
   dx, dy = mesh.dx, mesh.dy

   for i=1:nx, j=1:ny
      dex_dy  = (fdtd.ex[i,mod1(j+1,ny)]-fdtd.ex[i,j]) / dy
      dey_dx  = (fdtd.ey[mod1(i+1,nx),j]-fdtd.ey[i,j]) / dx
      fdtd.bz[i,j] += dt * (dex_dy - dey_dx)
   end

   for i=1:nx, j=1:ny
      bz[i,j] = ( fdtd.bz[mod1(i-1,nx), mod1(j-1,ny)] + fdtd.bz[i, mod1(j-1,ny)] 
                + fdtd.bz[mod1(i-1,nx), j] + fdtd.bz[i, j]) / 4
   end
   
end 

function ampere_maxwell!( ex, ey, fdtd :: FDTD, mesh, dt)

   nx, ny = mesh.nx, mesh.ny
   dx, dy = mesh.dx, mesh.dy

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



function compute_current!( jx, jy, fdtd :: FDTD, p, mesh )


    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    
    fill!(jx, 0)
    fill!(jy, 0)
    
    for ipart=1:p.nbpart
    
       i = p.cell[ipart,1]
       j = p.cell[ipart,2]
    
       ip1 = mod1(i+1, nx)
       jp1 = mod1(j+1, ny)
    
       xp = p.pos[ipart,1]
       yp = p.pos[ipart,2]
    
       w = (mesh.dimx * mesh.dimy) / p.nbpart / (dx * dy)
    
       a1 = (mesh.x[i+1]-xp) * (mesh.y[j+1]-yp) * w
       a2 = (xp-mesh.x[i]) * (mesh.y[j+1]-yp) * w
       a3 = (xp-mesh.x[i]) * (yp-mesh.y[j]) * w
       a4 = (mesh.x[i+1]-xp) * (yp-mesh.y[j]) * w
    
       w = p.vit[ipart,1] / (dx*dy) 
    
       jx[i,j]     +=  a1*w  
       jx[ip1,j]   +=  a2*w 
       jx[ip1,jp1] +=  a3*w 
       jx[i,jp1]   +=  a4*w 
    
       w = p.vit[ipart,2] / (dx*dy) 
    
       jy[i,j]     += a1*w  
       jy[ip1,j]   += a2*w 
       jy[ip1,jp1] += a3*w 
       jy[i,jp1]   += a4*w 
    
    end
    
    for i=1:nx, j=1:ny
       fdtd.jx[i,j] = 0.5 * (jx[i,j]+jx[mod1(i+1,nx),j])
    end
    
    for i=1:nx, j=1:ny
       fdtd.jy[i,j] = 0.5 * (jy[i,j]+jy[i,mod1(j+1,ny)])
    end

end 

