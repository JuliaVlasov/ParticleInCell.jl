function calcul_j_cic( jx, jy, m, p, f )

   dx, dy = m.dx, m.dy

   fill!(jx, 0)
   fill!(jy, 0)
   
   for ipart=1:nbpart

      i = p.cell[ipart,1]
      j = p.cell[ipart,2]

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      dum = p.p[ipart] / (dx*dy)

      a1 = (x[i+1]-xp) * (y(j+1)-yp) * dum
      a2 = (xp-x[i]) * (y(j+1)-yp) * dum
      a3 = (xp-x[i]) * (yp-y[j]) * dum
      a4 = (x[i+1]-xp) * (yp-y[j]) * dum

      dum = p.vit[ipart,1] / (dx*dy) 

      jx[i,j]     = jx[i,j]     + a1*dum  
      jx[i+1,j]   = jx[i+1,j]   + a2*dum 
      jx[i+1,j+1] = jx[i+1,j+1] + a3*dum 
      jx[i,j+1]   = jx[i,j+1]   + a4*dum 

      dum = p.vit[ipart,2] / (dx*dy) 

      jy[i,j]     = jy[i,j]     + a1*dum  
      jy[i+1,j]   = jy[i+1,j]   + a2*dum 
      jy[i+1,j+1] = jy[i+1,j+1] + a3*dum 
      jy[i,j+1]   = jy[i,j+1]   + a4*dum 

   end

   for i=0:nx

      jx[i,0]  = jx[i,0] + jx[i,ny]
      jx[i,ny] = jx[i,0]
      jy[i,0]  = jy[i,0] + jy[i,ny]
      jy[i,ny] = jy[i,0]

   end

   for j=0:ny

      jx[0,j]  = jx[0,j] + jx[nx,j]
      jx[nx,j] = jx[0,j]
      jy[0,j]  = jy[0,j] + jy[nx,j]
      jy[nx,j] = jy[0,j]

   end


    for i=0:nx-1
    for j=0:ny
       f.jx[i,j] = 0.5 * (jx[i,j]+jx[i+1,j])
    end 
    end
    
    for i=0:nx
    for j=0:ny-1
       f.jy[i,j] = 0.5 * (jy[i,j]+jy[i,j+1])
    end
    end

end 
