export interpol_eb_m6!

function interpol_eb_m6!( e      :: Array{Float64,3}, 
                          fields :: MeshFields, 
                          x      :: Array{ComplexF64,3},
                          nbpart :: Int64,
                          ntau   :: Int64)

    nx, ny = fields.mesh.nx, fields.mesh.ny
    dx, dy = fields.mesh.dx, fields.mesh.dy

    xmin, xmax = fields.mesh.xmin, fields.mesh.xmax
    ymin, ymax = fields.mesh.ymin, fields.mesh.ymax
    dimx = xmax - xmin
    dimy = ymax - ymin

    for k=1:nbpart, n=1:ntau
    
       xn = mod( real(x[n,1,k]) - xmin, dimx)
       yn = mod( real(x[n,2,k]) - ymin, dimy)

       px = real(xn)/dx
       py = real(yn)/dy

       i = trunc(Int32, px)
       j = trunc(Int32, py)

       dpx = px - i
       dpy = py - j
    
       im3 = mod(i-3,nx) + 1
       im2 = mod(i-2,nx) + 1
       im1 = mod(i-1,nx) + 1
       ip1 = mod(i+1,nx) + 1
       ip2 = mod(i+2,nx) + 1
       ip3 = mod(i+3,nx) + 1
       jm3 = mod(j-3,ny) + 1
       jm2 = mod(j-2,ny) + 1
       jm1 = mod(j-1,ny) + 1
       jp1 = mod(j+1,ny) + 1
       jp2 = mod(j+2,ny) + 1
       jp3 = mod(j+3,ny) + 1

       i = i + 1
       j = j + 1

       cm3x = f_m6(3+dpx)
       cp3x = f_m6(3-dpx)
       cm2x = f_m6(2+dpx)
       cp2x = f_m6(2-dpx)
       cm1x = f_m6(1+dpx)
       cp1x = f_m6(1-dpx)
       cx   = f_m6(  dpx)
       cy   = f_m6(  dpy)
       cm3y = f_m6(3+dpy)
       cp3y = f_m6(3-dpy)
       cm2y = f_m6(2+dpy)
       cp2y = f_m6(2-dpy)
       cm1y = f_m6(1+dpy)
       cp1y = f_m6(1-dpy)
    
     
       for l = 1:2

           s = 0.0
           s += cm3x * cm3y * fields.e[l,im3,jm3]   
           s += cm3x * cm2y * fields.e[l,im3,jm2]   
           s += cm3x * cm1y * fields.e[l,im3,jm1]   
           s += cm3x * cy   * fields.e[l,im3,j  ]   
           s += cm3x * cp1y * fields.e[l,im3,jp1]   
           s += cm3x * cp2y * fields.e[l,im3,jp2]   
           s += cm3x * cp3y * fields.e[l,im3,jp3]   
           s += cm2x * cm3y * fields.e[l,im2,jm3]   
           s += cm2x * cm2y * fields.e[l,im2,jm2]   
           s += cm2x * cm1y * fields.e[l,im2,jm1]   
           s += cm2x * cy   * fields.e[l,im2,j  ]   
           s += cm2x * cp1y * fields.e[l,im2,jp1]   
           s += cm2x * cp2y * fields.e[l,im2,jp2]   
           s += cm2x * cp3y * fields.e[l,im2,jp3]   
           s += cm1x * cm3y * fields.e[l,im1,jm3]   
           s += cm1x * cm2y * fields.e[l,im1,jm2]   
           s += cm1x * cm1y * fields.e[l,im1,jm1]   
           s += cm1x * cy   * fields.e[l,im1,j  ]   
           s += cm1x * cp1y * fields.e[l,im1,jp1]   
           s += cm1x * cp2y * fields.e[l,im1,jp2]   
           s += cm1x * cp3y * fields.e[l,im1,jp3]   
           s += cx   * cm3y * fields.e[l,i  ,jm3]   
           s += cx   * cm2y * fields.e[l,i  ,jm2]   
           s += cx   * cm1y * fields.e[l,i  ,jm1]   
           s += cx   * cy   * fields.e[l,i  ,j  ]   
           s += cx   * cp1y * fields.e[l,i  ,jp1]   
           s += cx   * cp2y * fields.e[l,i  ,jp2]   
           s += cx   * cp3y * fields.e[l,i  ,jp3]   
           s += cp1x * cm3y * fields.e[l,ip1,jm3]   
           s += cp1x * cm2y * fields.e[l,ip1,jm2]   
           s += cp1x * cm1y * fields.e[l,ip1,jm1]   
           s += cp1x * cy   * fields.e[l,ip1,j  ]   
           s += cp1x * cp1y * fields.e[l,ip1,jp1]   
           s += cp1x * cp2y * fields.e[l,ip1,jp2]   
           s += cp1x * cp3y * fields.e[l,ip1,jp3]   
           s += cp2x * cm3y * fields.e[l,ip2,jm3]   
           s += cp2x * cm2y * fields.e[l,ip2,jm2]   
           s += cp2x * cm1y * fields.e[l,ip2,jm1]   
           s += cp2x * cy   * fields.e[l,ip2,j  ]   
           s += cp2x * cp1y * fields.e[l,ip2,jp1]   
           s += cp2x * cp2y * fields.e[l,ip2,jp2]   
           s += cp2x * cp3y * fields.e[l,ip2,jp3]   
           s += cp3x * cm3y * fields.e[l,ip3,jm3]   
           s += cp3x * cm2y * fields.e[l,ip3,jm2]   
           s += cp3x * cm1y * fields.e[l,ip3,jm1]   
           s += cp3x * cy   * fields.e[l,ip3,j  ]   
           s += cp3x * cp1y * fields.e[l,ip3,jp1]   
           s += cp3x * cp2y * fields.e[l,ip3,jp2]   
           s += cp3x * cp3y * fields.e[l,ip3,jp3]

           e[n,l,k] = s
    
       end

    end
    

end 

function interpol_eb_m6!( particles :: Particles, fields :: MeshFields )

    nx = fields.mesh.nx
    ny = fields.mesh.ny

    dx = fields.mesh.dx
    dy = fields.mesh.dy

    xmin, xmax = fields.mesh.xmin, fields.mesh.xmax
    ymin, ymax = fields.mesh.ymin, fields.mesh.ymax
    dimx = xmax - xmin
    dimy = ymax - ymin

    for k=1:particles.nbpart
    
       xn = mod( particles.x[1,k] - xmin, dimx)
       yn = mod( particles.x[2,k] - ymin, dimy)

       px = xn/dx
       py = yn/dy

       i = trunc(Int32, px)
       j = trunc(Int32, py)

       dpx = px - i
       dpy = py - j

       particles.x[1,k] = xmin + xn
       particles.x[2,k] = ymin + yn

       im3 = mod(i-3,nx) + 1
       im2 = mod(i-2,nx) + 1
       im1 = mod(i-1,nx) + 1
       ip1 = mod(i+1,nx) + 1
       ip2 = mod(i+2,nx) + 1
       ip3 = mod(i+3,nx) + 1
       jm3 = mod(j-3,ny) + 1
       jm2 = mod(j-2,ny) + 1
       jm1 = mod(j-1,ny) + 1
       jp1 = mod(j+1,ny) + 1
       jp2 = mod(j+2,ny) + 1
       jp3 = mod(j+3,ny) + 1

       i = i + 1
       j = j + 1

       cm3x = f_m6(3+dpx)
       cp3x = f_m6(3-dpx)
       cm2x = f_m6(2+dpx)
       cp2x = f_m6(2-dpx)
       cm1x = f_m6(1+dpx)
       cp1x = f_m6(1-dpx)
       cx   = f_m6(  dpx)
       cy   = f_m6(  dpy)
       cm3y = f_m6(3+dpy)
       cp3y = f_m6(3-dpy)
       cm2y = f_m6(2+dpy)
       cp2y = f_m6(2-dpy)
       cm1y = f_m6(1+dpy)
       cp1y = f_m6(1-dpy)
    
     
       for l = 1:2

           e = 0.0
           e += cm3x * cm3y * fields.e[l,im3,jm3]   
           e += cm3x * cm2y * fields.e[l,im3,jm2]   
           e += cm3x * cm1y * fields.e[l,im3,jm1]   
           e += cm3x * cy   * fields.e[l,im3,j  ]   
           e += cm3x * cp1y * fields.e[l,im3,jp1]   
           e += cm3x * cp2y * fields.e[l,im3,jp2]   
           e += cm3x * cp3y * fields.e[l,im3,jp3]   
           e += cm2x * cm3y * fields.e[l,im2,jm3]   
           e += cm2x * cm2y * fields.e[l,im2,jm2]   
           e += cm2x * cm1y * fields.e[l,im2,jm1]   
           e += cm2x * cy   * fields.e[l,im2,j  ]   
           e += cm2x * cp1y * fields.e[l,im2,jp1]   
           e += cm2x * cp2y * fields.e[l,im2,jp2]   
           e += cm2x * cp3y * fields.e[l,im2,jp3]   
           e += cm1x * cm3y * fields.e[l,im1,jm3]   
           e += cm1x * cm2y * fields.e[l,im1,jm2]   
           e += cm1x * cm1y * fields.e[l,im1,jm1]   
           e += cm1x * cy   * fields.e[l,im1,j  ]   
           e += cm1x * cp1y * fields.e[l,im1,jp1]   
           e += cm1x * cp2y * fields.e[l,im1,jp2]   
           e += cm1x * cp3y * fields.e[l,im1,jp3]   
           e += cx   * cm3y * fields.e[l,i  ,jm3]   
           e += cx   * cm2y * fields.e[l,i  ,jm2]   
           e += cx   * cm1y * fields.e[l,i  ,jm1]   
           e += cx   * cy   * fields.e[l,i  ,j  ]   
           e += cx   * cp1y * fields.e[l,i  ,jp1]   
           e += cx   * cp2y * fields.e[l,i  ,jp2]   
           e += cx   * cp3y * fields.e[l,i  ,jp3]   
           e += cp1x * cm3y * fields.e[l,ip1,jm3]   
           e += cp1x * cm2y * fields.e[l,ip1,jm2]   
           e += cp1x * cm1y * fields.e[l,ip1,jm1]   
           e += cp1x * cy   * fields.e[l,ip1,j  ]   
           e += cp1x * cp1y * fields.e[l,ip1,jp1]   
           e += cp1x * cp2y * fields.e[l,ip1,jp2]   
           e += cp1x * cp3y * fields.e[l,ip1,jp3]   
           e += cp2x * cm3y * fields.e[l,ip2,jm3]   
           e += cp2x * cm2y * fields.e[l,ip2,jm2]   
           e += cp2x * cm1y * fields.e[l,ip2,jm1]   
           e += cp2x * cy   * fields.e[l,ip2,j  ]   
           e += cp2x * cp1y * fields.e[l,ip2,jp1]   
           e += cp2x * cp2y * fields.e[l,ip2,jp2]   
           e += cp2x * cp3y * fields.e[l,ip2,jp3]   
           e += cp3x * cm3y * fields.e[l,ip3,jm3]   
           e += cp3x * cm2y * fields.e[l,ip3,jm2]   
           e += cp3x * cm1y * fields.e[l,ip3,jm1]   
           e += cp3x * cy   * fields.e[l,ip3,j  ]   
           e += cp3x * cp1y * fields.e[l,ip3,jp1]   
           e += cp3x * cp2y * fields.e[l,ip3,jp2]   
           e += cp3x * cp3y * fields.e[l,ip3,jp3]

           particles.e[l,k] = e
    
       end

    end
    

end 
