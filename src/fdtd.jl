export FDTD

struct FDTD

    m :: Mesh
    ebj :: Array{Float64, 3}
    ex :: Array{Float64,2}
    ey :: Array{Float64,2}
    bz :: Array{Float64,2}
    jx :: Array{Float64,2}
    jy :: Array{Float64,2}

    function FDTD( mesh )

        nx, ny = mesh.nx, mesh.ny
        ebj = zeros(5, nx+1,ny+1)
        ex = zeros(nx,ny+1)
        ey = zeros(nx+1,ny)
        bz = zeros(nx,ny)
        jx = zeros(nx,ny+1)
        jy = zeros(nx+1,ny)

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

   for i=1:nx+1, j=1:ny+1
      fdtd.ebj[3,i,j] = ( fdtd.bz[mod1(i-1,nx), mod1(j-1,ny)] 
                        + fdtd.bz[mod1(i,  nx), mod1(j-1,ny)] 
                        + fdtd.bz[mod1(i-1,nx), mod1(j,  ny)] 
                        + fdtd.bz[mod1(i,  nx), mod1(j,  ny)]) / 4
   end

end 

export ampere_maxwell!

function ampere_maxwell!( fdtd :: FDTD, dt)

   nx, ny = fdtd.m.nx, fdtd.m.ny
   dx, dy = fdtd.m.dx, fdtd.m.dy

   for i=1:nx, j=1:ny+1
      dbz_dy = (fdtd.bz[i,mod1(j,ny)]-fdtd.bz[i,mod1(j-1,ny)]) / dy
      fdtd.ex[i,j] += dt * dbz_dy - dt * fdtd.jx[i,j]
   end

   for i=1:nx+1, j=1:ny
      dbz_dx = (fdtd.bz[mod1(i,nx),mod1(j,ny)]-fdtd.bz[mod1(i-1,nx),j]) / dx
      fdtd.ey[i,j] -= dt * dbz_dx - dt * fdtd.jy[i,j]
   end

   for i=1:nx+1, j=1:ny+1
      fdtd.ebj[1,i,j] = 0.5*( fdtd.ex[mod1(i-1,nx),j] + fdtd.ex[mod1(i,nx),j] )
      fdtd.ebj[2,i,j] = 0.5*( fdtd.ey[i,mod1(j-1,ny)] + fdtd.ey[i,mod1(j,ny)])
   end

end 

export interpol_eb!

function interpol_eb!(p :: Array{Float64,2}, fdtd::FDTD)

    nbpart = size(p)[2]
    dx = fdtd.m.dx
    dy = fdtd.m.dy

    @inbounds for ipart = 1:nbpart

        xp = p[1, ipart] / dx
        yp = p[2, ipart] / dy

        i = floor(Int, xp ) + 1
        j = floor(Int, yp ) + 1

        dxp = xp - i + 1
        dyp = yp - j + 1
    
        a1 = (1-dxp) * (1-dyp)
        a2 = dxp * (1-dyp)
        a3 = dxp * dyp
        a4 = (1-dxp) * dyp

        p[5,ipart] = a1 * fdtd.ebj[1, i, j] + a2 * fdtd.ebj[1, i+1, j] + a3 * fdtd.ebj[1, i+1, j+1] + a4 * fdtd.ebj[1, i, j+1]
        p[6,ipart] = a1 * fdtd.ebj[2, i, j] + a2 * fdtd.ebj[2, i+1, j] + a3 * fdtd.ebj[2, i+1, j+1] + a4 * fdtd.ebj[2, i, j+1]
        p[7,ipart] = a1 * fdtd.ebj[3, i, j] + a2 * fdtd.ebj[3, i+1, j] + a3 * fdtd.ebj[3, i+1, j+1] + a4 * fdtd.ebj[3, i, j+1]

    end


end

export compute_current!

function compute_current!( fdtd :: FDTD, p )

    nbpart = size(p)[2]
    mesh = fdtd.m
    nx, ny = fdtd.m.nx, fdtd.m.ny
    dx, dy = fdtd.m.dx, fdtd.m.dy
    
    fdtd.ebj[4,:,:] .= 0
    fdtd.ebj[5,:,:] .= 0

    factor = ( nx * ny ) / nbpart
    
    for ipart=1:nbpart
    
       xp = p[1,ipart] / dx
       yp = p[2,ipart] / dy

       i = floor(Int, xp) + 1
       j = floor(Int, yp) + 1
    
       dxp = xp - i + 1
       dyp = yp - j + 1
    
       a1 = (1-dxp) * (1-dyp)
       a2 = dxp * (1-dyp)
       a3 = dxp * dyp
       a4 = (1-dxp) * dyp

       w1 = p[3,ipart] * factor
       w2 = p[4,ipart] * factor
    
       fdtd.ebj[4,i,j]     += a1*w1
       fdtd.ebj[5,i,j]     += a1*w2  

       fdtd.ebj[4,i+1,j]   += a2*w1 
       fdtd.ebj[5,i+1,j]   += a2*w2 

       fdtd.ebj[4,i+1,j+1] += a3*w1 
       fdtd.ebj[5,i+1,j+1] += a3*w2 

       fdtd.ebj[4,i,j+1]   += a4*w1 
       fdtd.ebj[5,i,j+1]   += a4*w2 
    
    end

    for i=1:nx+1
      fdtd.ebj[4:5,i,1]  .+= fdtd.ebj[4:5,i,ny+1]
      fdtd.ebj[4:5,i,ny+1]  .= fdtd.ebj[4:5,i,1]
    end
    for j=1:ny+1
      fdtd.ebj[4:5,1,j]  .+= fdtd.ebj[4:5,nx+1,j]
      fdtd.ebj[4:5,nx+1,j]  .= fdtd.ebj[4:5,1,j]
    end

    for i=1:nx, j=1:ny+1
       fdtd.jx[i,j] = 0.5 * (fdtd.ebj[4,i,j]+fdtd.ebj[4,i+1,j])
    end
    
    for i=1:nx+1, j=1:ny
       fdtd.jy[i,j] = 0.5 * (fdtd.ebj[5,i,j]+fdtd.ebj[5,i,j+1])
    end

end 

export compute_energy

compute_energy( fdtd :: FDTD ) = 0.5 * log( sum(view(fdtd.ebj,1,:,:).^2) * fdtd.m.dx * fdtd.m.dy)
