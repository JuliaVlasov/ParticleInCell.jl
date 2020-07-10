"""
Calcul des composantes des champs 
sur les noeuds du maillage de rho
par interpolation lineaire
"""
function decalage!( m :: Mesh, f :: Fields, ex, ey, bz )

   nx, ny = f.nx, f.ny

    for i=1:nx-1
       for j=1:ny-1
          ex[i,j] = ( m.hx[i]*f.ex[i-1,j] + m.hx[i-1]*f.ex[i,j] ) / (m.hx[i]+m.hx[i-1])
          ey[i,j] = ( m.hy[j]*f.ey[i,j-1] + m.hy[j-1]*f.ey[i,j] ) / (m.hy[j]+m.hy[j-1])
          bz[i,j] = ((( m.hx[i]*f.bz[i-1,j-1] + m.hx[i-1]*f.bz[i,j-1] )
                * m.hy[j] + ( m.hx[i]*f.bz[i-1,j] + m.hx[i-1]*f.bz[i,j] ) 
                * m.hy[j-1] ) / ( (m.hx[i]+m.hx[i-1]) * (m.hy[j]+m.hy[j-1]) ))
       end
    end

   for i=1:nx-1
      ex[i,0] = ( m.hx[i]*f.ex[i-1,0] + m.hx[i-1]*f.ex[i,0] ) / (m.hx[i]+m.hx[i-1])
      ey[i,0] = ( m.hy[0]*f.ey[i,ny-1] + m.hy[ny-1]*f.ey[i,0] ) / (m.hy[0]+m.hy[ny-1])
      bz[i,0] = ( ( m.hx[i]*f.bz[i-1,ny-1] + m.hx[i-1]*f.bz[i,ny-1] )* m.hy[0] + ( m.hx[i]*f.bz[i-1,0] + m.hx[i-1]*f.bz[i,0] ) * m.hy[ny-1] ) / ( (m.hx[i]+m.hx[i-1]) * (m.hy[0]+m.hy[ny-1]) )
      ex[i,ny] = ex[i,0] 
      ey[i,ny] = ey[i,0] 
      bz[i,ny] = bz[i,0] 
   end

   for j=1:ny-1 
      ex[0,j] = ( m.hx[0]*f.ex[nx-1,j] + m.hx[nx-1]*f.ex[0,j] ) / (m.hx[0]+m.hx[nx-1])
      ey[0,j] = ( m.hy[j]*f.ey[0,j-1] + m.hy[j-1]*f.ey[0,j] ) / (m.hy[j]+m.hy[j-1])
      bz[0,j] = ( ( m.hx[0]*f.bz[nx-1,j-1] + m.hx[nx-1]*f.bz[0,j-1] ) * m.hy[j] + ( m.hx[0]*f.bz[nx-1,j] + m.hx[nx-1]*f.bz[0,j] ) * m.hy[j-1] ) / ( (m.hx[0]+m.hx[nx-1]) * (m.hy[j]+m.hy[j-1]) )
      ex[nx,j] = ex[0,j] 
      ey[nx,j] = ey[0,j] 
      bz[nx,j] = bz[0,j] 
   end

   ex[0,0] = ( m.hx[0]*f.ex[nx-1,0] + m.hx[nx-1]*f.ex[0,0] ) / (m.hx[0]+m.hx[nx-1])
   ey[0,0] = ( m.hy[0]*f.ey[0,ny-1] + m.hy[ny-1]*f.ey[0,0] ) / (m.hy[0]+m.hy[ny-1])
   bz[0,0] = ( ( m.hx[0]*f.bz[nx-1,ny-1] + m.hx[nx-1]*f.bz[0,ny-1] ) * m.hy[0] + ( m.hx[0]*f.bz[nx-1,0] + m.hx[nx-1]*f.bz[0,0] ) * m.hy[ny-1] ) / ( (m.hx[0]+m.hx[nx-1]) * (m.hy[0]+m.hy[ny-1]) )

   ex[nx,0]  = ex[0,0] 
   ex[nx,ny] = ex[0,0] 
   ex[0,ny]  = ex[0,0] 

   ey[nx,0]  = ey[0,0] 
   ey[nx,ny] = ey[0,0] 
   ey[0,ny]  = ey[0,0] 

   bz[nx,0]  = bz[0,0] 
   bz[nx,ny] = bz[0,0] 
   bz[0,ny]  = bz[0,0] 

end
