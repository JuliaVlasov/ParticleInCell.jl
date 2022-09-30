function compute_rho_cic!( rho, mesh, x, w)

dx = mesh.dx
dy = mesh.dy
dz = mesh.dz

fill!(rho, 0.0)

vol = w / (dx*dy*dz)

nbpart = size(x)[2]

for m in 1:nbpart

   xp = x[1, m]/dx
   yp = x[2, m]/dy
   zp = x[3, m]/dz

   ip = floor(Int, xp)
   jp = floor(Int, yp)
   kP = floor(Int, zp)

   dxp = xp - ip
   dyp = yp - jp
   dzp = zp - kp

   a1 = (1 - dxp) * (1 - dyp) * (1 - dzp) 
   a2 = dxp * (1 - dyp) * (1 - dzp) 
   a3 = (1 - dxp) * dyp * (1 - dzp) 
   a4 = dxp * dyp * (1 - dzp)
   a5 = (1 - dxp) * (1 - dyp) * dzp
   a6 = dxp * (1 - dyp) * dzp
   a7 = (1 - dxp) * dyp * dzp
   a8 = dxp * dyp * dzp

   ip = ip+1
   jp = jp+1
   kp = kp+1

   rho[ip  ,jp  ,kp  ] += a1 * vol
   rho[ip+1,jp  ,kp  ] += a2 * vol
   rho[ip  ,jp+1,kp  ] += a3 * vol
   rho[ip+1,jp+1,kp  ] += a4 * vol
   rho[ip  ,jp  ,kp+1] += a5 * vol
   rho[ip+1,jp  ,kp+1] += a6 * vol
   rho[ip  ,jp+1,kp+1] += a7 * vol
   rho[ip+1,jp+1,kp+1] += a8 * vol

end

nx = mesh.nx
ny = mesh.ny
nz = mesh.nz

rho[nx+1,:,:] = rho[1,:,:]
rho[:,ny+1,:] = rho[:,1,:]
rho[:,:,nz+1] = rho[:,:,1]


end 
