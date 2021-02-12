export FDTD

struct FDTD

    m :: Mesh
    ex
    ey
    bz
    jx
    jy
    ebj
    

    function FDTD( m )

        nx, ny = m.nx, m.ny
        ex = zeros(nx,ny)
        ey = zeros(nx,ny)
        bz = zeros(nx,ny)
        jx = zeros(nx,ny)
        jy = zeros(nx,ny)
        ebj = zeros(3, nx, ny)

        new( m, ex, ey, bz, jx, jy, ebj )

    end

end

function faraday!( fdtd :: FDTD, dt )

   nx, ny = fdtd.m.nx, fdtd.m.ny
   dx, dy = fdtd.m.dx, fdtd.m.dy

   for j=1:ny, i=1:nx
      dex_dy  = (fdtd.ex[i,mod1(j+1,ny)]-fdtd.ex[i,j]) / dy
      dey_dx  = (fdtd.ey[mod1(i+1,nx),j]-fdtd.ey[i,j]) / dx
      fdtd.bz[i,j] += dt * (dex_dy - dey_dx)
   end

   for j=1:ny, i=1:nx
      fdtd.ebj[3,i,j] = ( fdtd.bz[mod1(i-1,nx), mod1(j-1,ny)] + fdtd.bz[i, mod1(j-1,ny)] 
                  + fdtd.bz[mod1(i-1,nx), j] + fdtd.bz[i, j]) / 4
   end
   
end 

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


function interpol_eb!(p::Particles, fdtd::FDTD)

    nx, ny = fdtd.m.nx, fdtd.m.ny
    dx, dy = fdtd.m.dx, fdtd.m.dy

    @inbounds for ip = 1:p.nbpart

        xp = p.data[1, ip] / dx
        yp = p.data[2, ip] / dy

        i = trunc(Int, xp ) + 1
        j = trunc(Int, yp ) + 1

        ip1 = mod1(i + 1, nx)
        jp1 = mod1(j + 1, ny)

        dxp = i - 1 - xp
        dyp = j - 1 - yp
        dxq = 1 - dxp
        dyq = 1 - dyp

        a1 = dxp * dyp 
        a2 = dxq * dyp
        a3 = dxq * dyq
        a4 = dxp * dyq

        p.data[5,ip] = a1 * fdtd.ebj[1,i,j] + a2 * fdtd.ebj[1,ip1,j] + a3 * fdtd.ebj[1,ip1,jp1] + a4 * fdtd.ebj[1,i,jp1]
        p.data[6,ip] = a1 * fdtd.ebj[2,i,j] + a2 * fdtd.ebj[2,ip1,j] + a3 * fdtd.ebj[2,ip1,jp1] + a4 * fdtd.ebj[2,i,jp1]
        p.data[7,ip] = a1 * fdtd.ebj[3,i,j] + a2 * fdtd.ebj[3,ip1,j] + a3 * fdtd.ebj[3,ip1,jp1] + a4 * fdtd.ebj[3,i,jp1]

    end

end

function compute_current!( fdtd :: FDTD, p :: Particles )


    nx, ny = fdtd.m.nx, fdtd.m.ny
    dx, dy = fdtd.m.dx, fdtd.m.dy
    
    fill!(fdtd.ebj, 0)
    
    for ipart=1:p.nbpart
    
       xp = p.data[1,ipart] / dx 
       yp = p.data[2,ipart] / dy

       i = trunc(Int, xp ) + 1
       j = trunc(Int, yp ) + 1
    
       ip1 = mod1(i+1, nx)
       jp1 = mod1(j+1, ny)
    
       dxp = i + 1 - xp
       dyp = j + 1 - yp
       dxq = 1 - dxp
       dyq = 1 - dyp

       a1 = dxp * dyp
       a2 = dxq * dyp
       a3 = dxq * dyq
       a4 = dxp * dyq
    
       w1 = p.data[3,ipart]
       w2 = p.data[4,ipart]
    
       fdtd.ebj[1,i,j]     += a1*w1  
       fdtd.ebj[2,i,j]     += a1*w2  

       fdtd.ebj[1,ip1,j]   += a2*w1
       fdtd.ebj[2,ip1,j]   += a2*w2

       fdtd.ebj[1,ip1,jp1] += a3*w1
       fdtd.ebj[2,ip1,jp1] += a3*w2

       fdtd.ebj[1,i,jp1]   += a4*w1
       fdtd.ebj[2,i,jp1]   += a4*w2
    
    end

    fdtd.ebj .*= (nx * ny) / p.nbpart 
    
    for i=1:nx, j=1:ny
       fdtd.jx[i,j] = 0.5 * (fdtd.ebj[1,i,j]+fdtd.ebj[1,mod1(i+1,nx),j])
    end
    
    for i=1:nx, j=1:ny
       fdtd.jy[i,j] = 0.5 * (fdtd.ebj[2,i,j]+fdtd.ebj[2,i,mod1(j+1,ny)])
    end

end 

export compute_energy

compute_energy( fdtd :: FDTD ) = 0.5 * log( sum(fdtd.ex.^2) * fdtd.m.dx * fdtd.m.dy)

