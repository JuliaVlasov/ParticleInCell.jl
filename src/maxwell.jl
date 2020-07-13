export MaxwellSolver


"""
     MaxwellSolver(mesh, c, e0)

Initialize a Maxwell solver that uses FDTD numerical method.

- `c` : speed of light
- `e0` :  electric permitivity

"""
struct MaxwellSolver

   m :: Mesh
   ex
   ey
   bz
   jx
   jy
   c
   e0

   function MaxwellSolver(m, c, e0)

       nx = m.nx
       ny = m.ny

       ex = OffsetArray(zeros(nx  ,ny+1), 0:nx-1,0:ny)
       ey = OffsetArray(zeros(nx+1,ny  ), 0:nx  ,0:ny-1)
       bz = OffsetArray(zeros(nx  ,ny  ), 0:nx-1,0:ny-1)
       jx = OffsetArray(zeros(nx  ,ny+1), 0:nx-1,0:ny)
       jy = OffsetArray(zeros(nx+1,ny  ), 0:nx  ,0:ny-1)

       new( m, ex, ey, bz, jx, jy, c, e0)

   end

end

export faraday!

"""
    faraday!( f , bz, dt )

Solve the Faraday equation to compute the magnetic field at ``t+1/2`` from electric field 
at ``t`` and magnetic field at ``t-1/2``.
"""
function faraday!( f :: MaxwellSolver, bz, dt )

   nx, ny = f.m.nx, f.m.ny

   for j=0:ny-1, i=0:nx-1
      dex_dy    = (f.ex[i,j+1] - f.ex[i,j]) / f.m.hy[j]
      dey_dx    = (f.ey[i+1,j] - f.ey[i,j]) / f.m.hx[i]
      f.bz[i,j] = f.bz[i,j] + dt * (dex_dy - dey_dx)
   end

   for j=1:ny-1
   for i=1:nx-1
       bz[i,j] = ((( f.m.hx[i]*f.bz[i-1,j-1] + f.m.hx[i-1]*f.bz[i,j-1] )
                * f.m.hy[j] + ( f.m.hx[i]*f.bz[i-1,j] + f.m.hx[i-1]*f.bz[i,j] ) 
                * f.m.hy[j-1] ) / ( (f.m.hx[i]+f.m.hx[i-1]) * (f.m.hy[j]+f.m.hy[j-1]) ))
   end
   end

   for i=1:nx-1
      bz[i,0] = ( ( f.m.hx[i]*f.bz[i-1,ny-1] + f.m.hx[i-1]*f.bz[i,ny-1] )* f.m.hy[0] + ( f.m.hx[i]*f.bz[i-1,0] + f.m.hx[i-1]*f.bz[i,0] ) * f.m.hy[ny-1] ) / ( (f.m.hx[i]+f.m.hx[i-1]) * (f.m.hy[0]+f.m.hy[ny-1]) )
      bz[i,ny] = bz[i,0] 
   end

   bz[0,0] = ( ( f.m.hx[0]*f.bz[nx-1,ny-1] + f.m.hx[nx-1]*f.bz[0,ny-1] ) * f.m.hy[0] + ( f.m.hx[0]*f.bz[nx-1,0] + f.m.hx[nx-1]*f.bz[0,0] ) * f.m.hy[ny-1] ) / ( (f.m.hx[0]+f.m.hx[nx-1]) * (f.m.hy[0]+f.m.hy[ny-1]) )
   bz[nx,0]  = bz[0,0] 
   bz[nx,ny] = bz[0,0] 
   bz[0,ny]  = bz[0,0] 


end 

export ampere_maxwell!

"""
     ampere_maxwell!( f :: MaxwellSolver, ex, ey, dt )

compute the electric field vector at t+1 from electrict field at t and magnetic field at
t-1/2. Solutions are compute on staggered grid (Ex[i+1/2,j] and Ey[i,j+1/2])
then interpolate. 
"""
function ampere_maxwell!( f :: MaxwellSolver, ex, ey, dt )

   @show csq = f.c * f.c
   nx, ny = f.m.nx, f.m.ny

   for j=1:ny-1
   for i=0:nx-1
      dbz_dy = (f.bz[i,j]-f.bz[i,j-1]) / f.m.hhy[j]
      f.ex[i,j] = f.ex[i,j] + csq * dt * dbz_dy - dt * f.jx[i,j]/f.e0
   end
   end

   for j=0:ny-1
   for i=1:nx-1
      dbz_dx = (f.bz[i,j]-f.bz[i-1,j]) / f.m.hhx[i]
      f.ey[i,j] = f.ey[i,j] - csq * dt * dbz_dx - dt * f.jy[i,j]/f.e0
   end
   end

   for i=0:nx-1
       dbz_dy = (f.bz[i,0]-f.bz[i,ny-1]) / f.m.hhy[0]
       f.ex[i,0]  = f.ex[i,0] + csq * dt * dbz_dy - dt * f.jx[i,0]/f.e0
       f.ex[i,ny] = f.ex[i,0]
   end
    
   for j=0:ny-1
       dbz_dx = (f.bz[0,j]-f.bz[nx-1,j]) / f.m.hhx[0]
       f.ey[0,j]  = f.ey[0,j] - csq * dt * dbz_dx - dt * f.jy[0,j]/f.e0
       f.ey[nx,j] = f.ey[0,j]
   end

   for j=1:ny-1
   for i=1:nx-1
       ex[i,j] = ( f.m.hx[i]*f.ex[i-1,j] + f.m.hx[i-1]*f.ex[i,j] ) / (f.m.hx[i]+f.m.hx[i-1])
       ey[i,j] = ( f.m.hy[j]*f.ey[i,j-1] + f.m.hy[j-1]*f.ey[i,j] ) / (f.m.hy[j]+f.m.hy[j-1])
   end
   end

   for i=1:nx-1
      ex[i,0] = ( f.m.hx[i]*f.ex[i-1,0] + f.m.hx[i-1]*f.ex[i,0] ) / (f.m.hx[i]+f.m.hx[i-1])
      ey[i,0] = ( f.m.hy[0]*f.ey[i,ny-1] + f.m.hy[ny-1]*f.ey[i,0] ) / (f.m.hy[0]+f.m.hy[ny-1])
      ex[i,ny] = ex[i,0] 
      ey[i,ny] = ey[i,0] 
   end

   ex[0,0] = ( f.m.hx[0]*f.ex[nx-1,0] + f.m.hx[nx-1]*f.ex[0,0] ) / (f.m.hx[0]+f.m.hx[nx-1])
   ey[0,0] = ( f.m.hy[0]*f.ey[0,ny-1] + f.m.hy[ny-1]*f.ey[0,0] ) / (f.m.hy[0]+f.m.hy[ny-1])

   ex[nx,0]  = ex[0,0] 
   ex[nx,ny] = ex[0,0] 
   ex[0,ny]  = ex[0,0] 

   ey[nx,0]  = ey[0,0] 
   ey[nx,ny] = ey[0,0] 
   ey[0,ny]  = ey[0,0] 


end

