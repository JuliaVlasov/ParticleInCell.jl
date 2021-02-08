using Statistics

export compute_rho

function compute_rho( p, m )

   nx, ny = m.nx, m.ny
   dx, dy = m.dx, m.dy
   rho = zeros(nx, ny)
   
   for ipart=1:p.nbpart

      i = p.cell[ipart,1]
      j = p.cell[ipart,2]

      ip1 = mod1(i+1, nx)
      jp1 = mod1(j+1, ny)

      xp = p.pos[ipart,1]
      yp = p.pos[ipart,2]

      a1 = (m.x[i+1]-xp) * (m.y[j+1]-yp) 
      a2 = (xp-m.x[i]) * (m.y[j+1]-yp) 
      a3 = (xp-m.x[i]) * (yp-m.y[j]) 
      a4 = (m.x[i+1]-xp) * (yp-m.y[j]) 

      rho[i,j]     +=  a1 
      rho[ip1,j]   +=  a2 
      rho[ip1,jp1] +=  a3 
      rho[i,jp1]   +=  a4 

   end

   rho = rho .* (nx*ny) ./ p.nbpart ./ (dx*dy) 

   return rho .- mean(rho)

end 