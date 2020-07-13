# -*- coding: utf-8 -*-
include("../src/mesh.jl")
include("../src/maxwell.jl")

using Plots
using OffsetArrays


dimx, dimy = 1, 1
nx, ny = 64, 64
e0 = 1
c = 1
md, nd = 2, 2  # number of modes along each direction
dt = 0.001
nstep = 1 ÷ dt

mesh = Mesh( dimx, nx, dimy, ny )
maxwell = MaxwellSolver( mesh, c, e0) 
omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)

ex = OffsetArray(zeros(nx+1,ny+1), 0:nx,0:ny)
ey = OffsetArray(zeros(nx+1,ny+1), 0:nx,0:ny)
bz = OffsetArray(zeros(nx+1,ny+1), 0:nx,0:ny);

# Ex and Ey are set at t = 0.0
# Bz is set at  t = -dt/2

# +
for j=0:ny, i=0:nx
    bz[i,j] = - cos(md*pi*mesh.x[i]) * cos(nd*pi*mesh.y[j]) * cos(omega*(-0.5*dt))
end

maxwell.ex .= 0
maxwell.ey .= 0

faraday!(maxwell, bz, dt) 
extrema(bz)
contourf(mesh.x[0:nx], mesh.y[0:ny], bz[0:nx,0:ny],  
    aspect_ratio=:equal) # set indices to create a copy

# +
time = 0
dt

ampere_maxwell!(maxwell, ex, ey, dt) 

extrema(ex)
# -



# +
for istep = 1:nstep # Loop over time
    
    

    time = time + 0.5dt
    
    

    time = time + 0.5dt

end # next time step
# -

err_l2 = 0.0
time = (nstep-0.5)*dt
for j = 0:ny, i = 0:nx
    th_bz = (- cos(md*pi*mesh.x[i]) * cos(nd*pi*mesh.y[j]) * cos(omega*time))
    err_l2 += (bz[i,j] - th_bz)^2
end


err_l2


