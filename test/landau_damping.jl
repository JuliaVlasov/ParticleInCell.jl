@testset "Landau Damping" begin

    nx        = 128  # nombre de pts suivant x
    ny        = 16   # nombre de pts suivant y
    cfl       = 0.9  # nombre de Courant-Friedrich-Levy
    tfinal    = 50   # temps final
    c         = 8    # vitesse de la lumiere
    e0        = 1    # permittivite du vide
    
    alpha = 0.1
    kx = 0.5
    ky = 0.
    dimx = 2*pi/kx
    dimy = 1  
    poids = dimx * dimy 
    
    mesh = Mesh( dimx, nx, dimy, ny)
    
    dx = mesh.dx
    dy = mesh.dy
    dt = cfl  / sqrt(1/(dx*dx)+1/(dy*dy)) / c
    
    println(" cfl = $cfl ")
    println(" dx = $dx dy = $dy dt = $dt")
    
    maxwell = Maxwell(mesh)
    
    ex = zeros(nx,ny)
    ey = zeros(nx,ny)
    bz = zeros(nx,ny)
    jx = zeros(nx,ny)
    jy = zeros(nx,ny)
    rho = zeros(nx, ny)
    
    time  = 0
    
    for i=1:nx
        aux1 = alpha/kx * sin(kx*mesh.x[i])
        aux2 = alpha * cos(kx*mesh.x[i])
        for j=1:ny
            ex[i,j] = aux1
            rho[i,j] = aux2
        end
    end
          
    nbpart = 100*nx*ny
    
    particles = Particles(nbpart)
    
    landau_sampling!( particles, alpha, kx )
    update_cells!( particles, mesh )
    
    time = 0
    for istep in 1:10
    
       if time > 0
           faraday!( bz, maxwell, ex, ey, 0.5dt ) 
       end

       interpol_eb!( ex, ey, bz, particles, mesh )
    
       push_v!( particles, dt )
       push_x!( particles, 0.5dt)  # x(n) --> x(n+1/2)
       update_cells!( particles, mesh )
       compute_current!( jx, jy, particles, mesh)
       push_x!( particles, 0.5dt)  # x(n+1/2) -- x(n+1)
       update_cells!( particles, mesh )
       faraday!(bz, maxwell, ex, ey, 0.5dt)
       ampere_maxwell!(ex, ey, maxwell, bz, jx, jy, dt)
    
       time = time + dt
    
    end 

    @test true

end
