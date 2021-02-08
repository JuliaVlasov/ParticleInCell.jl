function compute_current( jx, jy, p, m )

   dx, dy = m.dx, m.dy

   fill!(jx, 0)
   fill!(jy, 0)
   
   for ipart=1:nbpart


      i = p.cell[ipart,1]
      j = p.cell[ipart,2]

      ip1 = mod1(i+1,nx)
      jp1 = mod1(j+1,ny)

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      a1 = (x[i+1]-xp) * (y(j+1)-yp)
      a2 = (xp-x[i]) * (y(j+1)-yp) 
      a3 = (xp-x[i]) * (yp-y[j]) 
      a4 = (x[i+1]-xp) * (yp-y[j]) 

      dum = p.vit[ipart,1]

      jx[i,j]     = jx[i,j]     + a1*dum  
      jx[ip1,j]   = jx[ip1,j]   + a2*dum 
      jx[ip1,jp1] = jx[ip1,jp1] + a3*dum 
      jx[i,jp1]   = jx[i,jp1]   + a4*dum 

      dum = p.vit[ipart,2]

      jy[i,j]     = jy[i,j]     + a1*dum  
      jy[ip1,j]   = jy[ip1,j]   + a2*dum 
      jy[ip1,jp1] = jy[ip1,jp1] + a3*dum 
      jy[i,jp1]   = jy[i,jp1]   + a4*dum 

   end

   jx .= jx .* (nx*ny) ./ p.nbpart ./ (dx*dy) 
   jy .= jy .* (nx*ny) ./ p.nbpart ./ (dx*dy) 

end 
