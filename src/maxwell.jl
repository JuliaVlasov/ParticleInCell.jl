"""
   On utilise l'equation de Faraday sur un demi pas
   de temps pour le calcul du champ magnetique  Bz 
   a l'instant n puis n+1/2 apres deplacement des
   particules
"""
function faraday!( f :: Fields, m :: Mesh, dt )

   nx, ny = f.nx, f.ny

   for i=0:nx-1
   for j=0:ny-1
      dex_dy    = (f.ex[i,j+1] - f.ex[i,j]) / m.hy[j]
      dey_dx    = (f.ey[i+1,j] - f.ey[i,j]) / m.hx[i]
      f.bz[i,j] = f.bz[i,j] + 0.5 * dt * (dex_dy - dey_dx)
   end
   end

end 

"""
   Calcul du champ electrique E au temps n+1
   sur les points internes du maillage
   Ex aux points (i+1/2,j)
   Ey aux points (i,j+1/2)
"""
function ampere!( f :: Fields, m :: Mesh, dt, c, e0 )

   csq = c * c
   nx, ny = f.nx, f.ny

   for i=0:nx-1
   for j=1:ny-1
      dbz_dy = (f.bz[i,j]-f.bz[i,j-1]) / m.hhy[j]
      f.ex[i,j] = f.ex[i,j] + csq * dt * dbz_dy - dt * f.jx[i,j]/e0
   end
   end

   for i=1:nx-1
   for j=0:ny-1
      dbz_dx = (f.bz[i,j]-f.bz[i-1,j]) / m.hhx[i]
      f.ey[i,j] = f.ey[i,j] - csq * dt * dbz_dx - dt * f.jy[i,j]/e0
   end
   end

end

function conditions_limites( f :: Fields, m :: Mesh, dt, c, e0 )

    csq = c * c
    nx, ny = f.nx, f.ny

    for i=0:nx-1
       dbz_dy = (f.bz[i,0]-f.bz[i,ny-1]) / m.hhy[0]
       f.ex[i,0]  = f.ex[i,0] + csq * dt * dbz_dy - dt * f.jx[i,0]/e0
       f.ex[i,ny] = f.ex[i,0]
    end
    
    for j=0:ny-1
       dbz_dx = (f.bz[0,j]-f.bz[nx-1,j]) / m.hhx[0]
       f.ey[0,j]  = f.ey[0,j] - csq * dt * dbz_dx - dt * f.jy[0,j]/e0
       f.ey[nx,j] = f.ey[0,j]
    end

end 


