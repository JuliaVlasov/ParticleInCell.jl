using ParticleInCell

@testset "Particles" begin

    nx        = 128	    # nombre de pts suivant x
    ny        = 16   	# nombre de pts suivant y
    tfinal    = 50     	# temps final
    charge    = 1       # charge d'une macro particule
    masse     = 1       # masse d'une macro particule
    
    q_over_m = charge / masse
    
    alpha = 0.1
    kx = 0.5
    dimx = 2*pi/kx
    dimy = 1  
    
    mesh = Mesh( dimx, nx, dimy, ny)
    
    dx = mesh.dx
    dy = mesh.dy
    
    println(" dx = $dx dy = $dy ")
    
    ex  = zeros(nx+1,ny+1)
    rho = zeros(nx+1,ny+1)
    
    for i=1:nx+1
        aux1 = alpha/kx * sin(kx*mesh.x[i])
        aux2 = alpha * cos(kx*mesh.x[i])
        for j=1:ny+1
            ex[i,j] = aux1
            rho[i,j] = aux2
        end
    end

    @show sum(rho[1:nx,1:ny])

    nbpart = 100*nx*ny
    
    particles = Particles(nbpart)
    
    landau_sampling!( particles, alpha, kx )
    
    update_cells!( particles, mesh)

    @show sum(compute_rho(particles, mesh)[1:nx,1:ny] * dx * dy) 

    @test maximum(abs.(rho .- compute_rho(particles, mesh))) < 1e-6
    
    
end
