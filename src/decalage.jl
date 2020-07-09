"""
Calcul des composantes des champs 
sur les noeuds du maillage de rho
par interpolation lineaire
"""
function decalage!( f :: MeshFields, ex, ey, bz )

for i=1:nx-1
   for j=1:ny-1
      ex[i,j] = ( hx[i]*f.ex[i-1,j] + hx[i-1]*f.ex[i,j] ) / (hx[i]+hx[i-1])
      ey[i,j] = ( hy[j]*f.ey[i,j-1] + hy[j-1]*f.ey[i,j] ) / (hy[j]+hy[j-1])
      bz[i,j] = ((( hx[i]*f.bz[i-1,j-1] + hx[i-1]*f.bz[i,j-1] )
            * hy[j] + ( hx[i]*f.bz[i-1,j] + hx[i-1]*f.bz[i,j] ) 
            * hy[j-1] ) / ( (hx[i]+hx[i-1]) * (hy[j]+hy[j-1]) ))
   end
end

   for i=1:nx-1
      ex[i,0] = ( hx[i]*tm%ex[i-1,0] + hx[i-1]*tm%ex[i,0] ) / (hx[i]+hx[i-1])
      ey[i,0] = ( hy[0]*tm%ey[i,ny-1] + hy[ny-1]*tm%ey[i,0] ) / (hy[0]+hy[ny-1])
      bz[i,0] = ( ( hx[i]*tm%bz[i-1,ny-1] + hx[i-1]*tm%bz[i,ny-1] )* hy[0] + ( hx[i]*tm%bz[i-1,0] + hx[i-1]*tm%bz[i,0] ) * hy[ny-1] ) / ( (hx[i]+hx[i-1]) * (hy[0]+hy[ny-1]) )
      ex[i,ny] = ex[i,0] 
      ey[i,ny] = ey[i,0] 
      bz[i,ny] = bz[i,0] 
   end

   for j=1:ny-1 
      ex[0,j] = ( hx[0]*tm%ex[nx-1,j] + hx[nx-1]*tm%ex[0,j] ) / (hx[0]+hx[nx-1])
      ey[0,j] = ( hy[j]*tm%ey[0,j-1] + hy[j-1]*tm%ey[0,j] ) / (hy[j]+hy[j-1])
      bz[0,j] = ( ( hx[0]*tm%bz[nx-1,j-1] + hx[nx-1]*tm%bz[0,j-1] ) * hy[j] + ( hx[0]*tm%bz[nx-1,j] + hx[nx-1]*tm%bz[0,j] ) * hy[j-1] ) / ( (hx[0]+hx[nx-1]) * (hy[j]+hy[j-1]) )
      ex[nx,j] = ex[0,j] 
      ey[nx,j] = ey[0,j] 
      bz[nx,j] = bz[0,j] 
   end

   ex[0,0] = ( hx[0]*tm%ex[nx-1,0] + hx[nx-1]*tm%ex[0,0] ) / (hx[0]+hx[nx-1])
   ey[0,0] = ( hy[0]*tm%ey[0,ny-1] + hy[ny-1]*tm%ey[0,0] ) / (hy[0]+hy[ny-1])
   bz[0,0] = ( ( hx[0]*tm%bz[nx-1,ny-1] + hx[nx-1]*tm%bz[0,ny-1] ) * hy[0] + ( hx[0]*tm%bz[nx-1,0] + hx[nx-1]*tm%bz[0,0] ) * hy[ny-1] ) / ( (hx[0]+hx[nx-1]) * (hy[0]+hy[ny-1]) )

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
