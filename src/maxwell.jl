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

       ex = zeros(nx,ny+1)
       ey = zeros(nx+1,ny)
       bz = zeros(nx,ny)
       jx = zeros(nx,ny+1)
       jy = zeros(nx+1,ny)

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

   for j=1:ny, i=1:nx
      dex_dy    = (f.ex[i,j+1] - f.ex[i,j]) / f.m.dy
      dey_dx    = (f.ey[i+1,j] - f.ey[i,j]) / f.m.dx
      f.bz[i,j] = f.bz[i,j] + dt * (dex_dy - dey_dx)
   end

   for j=1:ny+1, i=1:nx+1
       bz[i,j] = 0.25 * ( f.bz[mod1(i-1,nx),mod1(j-1,ny)] 
                        + f.bz[mod1(i, nx),mod1(j-1,ny)] 
                        + f.bz[mod1(i-1,nx),mod1(j, ny)] 
                        + f.bz[mod1(i,nx),mod1(j,ny)] )
   end

end 

export ampere_maxwell!

"""
     ampere_maxwell!( f :: MaxwellSolver, ex, ey, dt )

compute the electric field vector at t+1 from electrict field at t and magnetic field at
t-1/2. Solutions are compute on staggered grid (Ex[i+1/2,j] and Ey[i,j+1/2])
then interpolate. 
"""
function ampere_maxwell!( f :: MaxwellSolver, ex, ey, dt )

   csq = f.c * f.c
   nx, ny = f.m.nx, f.m.ny

   for j=1:ny+1, i=1:nx
      dbz_dy = (f.bz[i,mod1(j,ny)]-f.bz[i,mod1(j-1,ny)]) / f.m.dy
      f.ex[i,j] = f.ex[i,j] + csq * dt * dbz_dy - dt * f.jx[i,j]/f.e0
   end

   for j=1:ny, i=1:nx+1
      dbz_dx = (f.bz[mod1(i,nx),j]-f.bz[mod1(i-1,nx),j]) / f.m.dx
      f.ey[i,j] = f.ey[i,j] - csq * dt * dbz_dx - dt * f.jy[i,j]/f.e0
   end

   for j=1:ny+1, i=1:nx+1
       ex[i,j] = 0.5 * ( f.ex[mod1(i-1,nx),j] + f.ex[mod1(i,nx),j] ) 
       ey[i,j] = 0.5 * ( f.ey[i,mod1(j-1,ny)] + f.ey[i,mod1(j,ny)] ) 
   end


end

