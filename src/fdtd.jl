export FDTD

struct FDTD

    m :: Mesh
    ebj :: Array{Float64, 3}
    ex 
    ey
    bz
    jx
    jy

    function FDTD( mesh )

        nx, ny = mesh.nx, mesh.ny
        ebj = zeros(5, nx,ny)
        ex = zeros(nx,ny)
        ey = zeros(nx,ny)
        bz = zeros(nx,ny)
        jx = zeros(nx,ny)
        jy = zeros(nx,ny)

        new( mesh, ebj, ex, ey, bz, jx, jy )

    end

end

export faraday!

function faraday!( fdtd :: FDTD, dt )

   nx, ny = fdtd.m.nx, fdtd.m.ny
   dx, dy = fdtd.m.dx, fdtd.m.dy

   for i=1:nx, j=1:ny
      dex_dy  = (fdtd.ex[i,mod1(j+1,ny)]-fdtd.ex[i,j]) / dy
      dey_dx  = (fdtd.ey[mod1(i+1,nx),j]-fdtd.ey[i,j]) / dx
      fdtd.bz[i,j] += dt * (dex_dy - dey_dx)
   end

   for i=1:nx, j=1:ny
      fdtd.ebj[3,i,j] = ( fdtd.bz[mod1(i-1,nx), mod1(j-1,ny)] + fdtd.bz[i, mod1(j-1,ny)] 
                + fdtd.bz[mod1(i-1,nx), j] + fdtd.bz[i, j]) / 4
   end
   
end 

export ampere_maxwell!

function ampere_maxwell!( fdtd :: FDTD, dt)

   nx, ny = fdtd.m.nx, fdtd.m.ny
   dx, dy = fdtd.m.dx, fdtd.m.dy

   for i=1:nx, j=1:ny
      dbz_dy = (fdtd.bz[i,j]-fdtd.bz[i,mod1(j-1,ny)]) / dy
      fdtd.ex[i,j] += dt * dbz_dy - dt * fdtd.jx[i,j]
   end

   for i=1:nx, j=1:ny
      dbz_dx = (fdtd.bz[i,j]-fdtd.bz[mod1(i-1,nx),j]) / dx
      fdtd.ey[i,j] -= dt * dbz_dx - dt * fdtd.jy[i,j]
   end

   for i=1:nx, j=1:ny
      fdtd.ebj[1,i,j] = 0.5*( fdtd.ex[mod1(i-1,nx),j] + fdtd.ex[i,j] )
      fdtd.ebj[2,i,j] = 0.5*( fdtd.ey[i,mod1(j-1,ny)] + fdtd.ey[i,j] )
   end

end 


export interpol_eb!

function interpol_eb!(p::Particles, fdtd::FDTD)

    nx, ny = fdtd.m.nx, fdtd.m.ny
    dx, dy = fdtd.m.dx, fdtd.m.dy
    dum = 1 / (dx * dy)

    @inbounds for ipart = 1:p.nbpart


        xp = p.data[1, ipart]
        yp = p.data[2, ipart]

        i = floor(Int, xp / dx) + 1
        j = floor(Int, yp / dy) + 1

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        a1 = (fdtd.m.x[i+1] - xp) * (fdtd.m.y[j+1] - yp) * dum
        a2 = (xp - fdtd.m.x[i]) * (fdtd.m.y[j+1] - yp) * dum
        a3 = (xp - fdtd.m.x[i]) * (yp - fdtd.m.y[j]) * dum
        a4 = (fdtd.m.x[i+1] - xp) * (yp - fdtd.m.y[j]) * dum

        p.data[5,ipart] = a1 * fdtd.ebj[1, i, j] + a2 * fdtd.ebj[1, ip1, j] + a3 * fdtd.ebj[1, ip1, jp1] + a4 * fdtd.ebj[1, i, jp1]
        p.data[6,ipart] = a1 * fdtd.ebj[2, i, j] + a2 * fdtd.ebj[2, ip1, j] + a3 * fdtd.ebj[2, ip1, jp1] + a4 * fdtd.ebj[2, i, jp1]
        p.data[7,ipart] = a1 * fdtd.ebj[3, i, j] + a2 * fdtd.ebj[3, ip1, j] + a3 * fdtd.ebj[3, ip1, jp1] + a4 * fdtd.ebj[3, i, jp1]

    end

end

export compute_current!

function compute_current!( fdtd :: FDTD, p )


    mesh = fdtd.m
    nx, ny = fdtd.m.nx, fdtd.m.ny
    dx, dy = fdtd.m.dx, fdtd.m.dy
    
    fdtd.ebj[4,:,:] .= 0
    fdtd.ebj[5,:,:] .= 0
    
    for ipart=1:p.nbpart
    
       xp = p.data[1,ipart]
       yp = p.data[2,ipart]

       i = floor(Int, xp / dx) + 1
       j = floor(Int, yp / dy) + 1
    
       ip1 = mod1(i+1, nx)
       jp1 = mod1(j+1, ny)
    
       w = (mesh.dimx * mesh.dimy) / p.nbpart / (dx * dy)
    
       a1 = (mesh.x[i+1]-xp) * (mesh.y[j+1]-yp) * w
       a2 = (xp-mesh.x[i]) * (mesh.y[j+1]-yp) * w
       a3 = (xp-mesh.x[i]) * (yp-mesh.y[j]) * w
       a4 = (mesh.x[i+1]-xp) * (yp-mesh.y[j]) * w
    
       w = p.data[3,ipart] / (dx*dy) 
    
       fdtd.ebj[4,i,j]     +=  a1*w  
       fdtd.ebj[4,ip1,j]   +=  a2*w 
       fdtd.ebj[4,ip1,jp1] +=  a3*w 
       fdtd.ebj[4,i,jp1]   +=  a4*w 
    
       w = p.data[4,ipart] / (dx*dy) 
    
       fdtd.ebj[5,i,j]     += a1*w  
       fdtd.ebj[5,ip1,j]   += a2*w 
       fdtd.ebj[5,ip1,jp1] += a3*w 
       fdtd.ebj[5,i,jp1]   += a4*w 
    
    end
    
    for i=1:nx, j=1:ny
       fdtd.jx[i,j] = 0.5 * (fdtd.ebj[4,i,j]+fdtd.ebj[4,mod1(i+1,nx),j])
    end
    
    for i=1:nx, j=1:ny
       fdtd.jy[i,j] = 0.5 * (fdtd.ebj[5,i,j]+fdtd.ebj[5,i,mod1(j+1,ny)])
    end

end 

export compute_energy

compute_energy( fdtd :: FDTD ) = 0.5 * log( sum(fdtd.ex.^2) * fdtd.m.dx * fdtd.m.dy)
