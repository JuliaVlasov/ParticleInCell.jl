export FDTD

struct FDTD

    ex
    ey
    bz
    jx
    jy

    function FDTD( mesh )

        nx, ny = mesh.nx, mesh.ny
        ex = zeros(nx-1,ny)
        ey = zeros(nx,ny-1)
        bz = zeros(nx-1,ny-1)
        jx = zeros(nx-1,ny)
        jy = zeros(nx,ny-1)

        new( ex, ey, bz, jx, jy )

    end

end

function faraday!( bz, fdtd :: FDTD, mesh, dt )

   nx, ny = mesh.nx, mesh.ny
   dx, dy = mesh.dx, mesh.dy

   for i=1:nx-1, j=1:ny-1
      dex_dy  = (fdtd.ex[i,j+1]-fdtd.ex[i,j]) / dy
      dey_dx  = (fdtd.ey[i+1,j]-fdtd.ey[i,j]) / dx
      fdtd.bz[i,j] += dt * (dex_dy - dey_dx)
   end

   for i=1:nx, j=1:ny
      bz[i,j] = 0.25 * ( fdtd.bz[mod1(i-1,nx-1), mod1(j-1,ny-1)] 
                       + fdtd.bz[mod1(i  ,nx-1), mod1(j-1,ny-1)] 
                       + fdtd.bz[mod1(i-1,nx-1), mod1(j  ,ny-1)] 
                       + fdtd.bz[mod1(i  ,nx-1), mod1(j  ,ny-1)])
   end
   
end 

function ampere_maxwell!( ex, ey, fdtd :: FDTD, mesh, dt)

   @show nx, ny = mesh.nx, mesh.ny
   @show dx, dy = mesh.dx, mesh.dy
   @show dt

   for i=1:nx-1, j=1:ny
      dbz_dy = (fdtd.bz[i,mod1(j,ny-1)]-fdtd.bz[i,mod1(j-1,ny-1)]) / dy
      fdtd.ex[i,j] += dt * dbz_dy #- dt * fdtd.jx[i,j]
   end

   for i=1:nx, j=1:ny-1
      dbz_dx = (fdtd.bz[mod1(i,nx-1),j]-fdtd.bz[mod1(i-1,nx-1),j]) / dx
      fdtd.ey[i,j] -= dt * dbz_dx #- dt * fdtd.jy[i,j]
   end

   for i=1:nx, j=1:ny
      im1 = mod1(i-1,nx-1)
      jm1 = mod1(j-1,ny-1)
      ex[i,j] = 0.5*( fdtd.ex[im1,j] + fdtd.ex[mod1(i,nx-1),j] )
      ey[i,j] = 0.5*( fdtd.ey[i,jm1] + fdtd.ey[i,mod1(j,ny-1)] )
   end

end 
