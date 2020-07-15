export compute_rho

function compute_rho( p, m )

   nx, ny = m.nx, m.ny
   dx, dy = m.dx, m.dy
   rho = zeros(nx+1, ny+1)
   
   for ipart=1:p.nbpart

      i = p.cell[ipart,1]
      j = p.cell[ipart,2]

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      dum = p.p[ipart] / ( dx * dy ) 

      a1 = (m.x[i+1]-xp) * (m.y[j+1]-yp) * dum
      a2 = (xp-m.x[i]) * (m.y[j+1]-yp) * dum
      a3 = (xp-m.x[i]) * (yp-m.y[j]) * dum
      a4 = (m.x[i+1]-xp) * (yp-m.y[j]) * dum

      rho[i,j]     +=  a1 
      rho[i+1,j]   +=  a2 
      rho[i+1,j+1] +=  a3 
      rho[i,j+1]   +=  a4 

   end


   rho[1:nx,ny+1] .+= rho[1:nx,1]
   rho[nx+1,1:ny] .+= rho[1,1:ny]
   rho[nx+1,ny+1]  += rho[1,1]
    
   rho[1:nx,1] .= rho[1:nx,ny+1]
   rho[1,1:ny] .= rho[nx+1,1:ny]
   rho[1,1]     = rho[nx+1,ny+1]

   rho ./= (dx*dy)

   rho_total  = sum(rho[1:nx,1:ny]) * dx * dy
    
   return rho .- rho_total / m.dimx / m.dimy

end 
