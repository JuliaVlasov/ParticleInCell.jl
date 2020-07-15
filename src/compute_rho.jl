export compute_rho

function compute_rho( p, m )

   nx, ny = m.nx, m.ny
   rho = zeros(nx+1, ny+1)
   
   for ipart=1:p.nbpart

      i = p.case[ipart,1]
      j = p.case[ipart,2]

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      dum = p.p[ipart] / (m.dx*m.dy)

      a1 = (m.x[i+1]-xp) * (m.y[j+1]-yp) * dum
      a2 = (xp-m.x[i]) * (m.y[j+1]-yp) * dum
      a3 = (xp-m.x[i]) * (yp-m.y[j]) * dum
      a4 = (m.x[i+1]-xp) * (yp-m.y[j]) * dum

      rho[i,j]     +=  a1 
      rho[i+1,j]   +=  a2 
      rho[i+1,j+1] +=  a3 
      rho[i,j+1]   +=  a4 

   end

   for i=1:nx+1
      rho[i,1]    += rho[i,ny+1]
      rho[i,ny+1]  = rho[i,1]
   end

   for j=2:ny
      rho[1,j]   += rho[nx+1,j]
      rho[nx+1,j] = rho[1,j]
   end

   return rho ./ (m.dx * m.dy)

end 
