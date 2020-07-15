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

    nbpart = 100*nx*ny
    
    particles = Particles(nbpart)
    
    landau_sampling!( particles, alpha, kx )
    
    xp = view(particles.pos, :, 1)
    vp = particles.vit
    
    # pp = plot(layout=(3,1))
    # histogram!(pp[1,1], xp, normalize=true, bins = 100, lab="x")
    # plot!(pp[1,1], x -> (1+alpha*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
    # histogram!(pp[2,1], vp[:,1], normalize=true, bins = 100, lab="vx")
    # plot!(pp[2,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
    # histogram!(pp[3,1], vp[:,2], normalize=true, bins = 100, lab="vy")
    # plot!(pp[3,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
    # savefig("particles.png")

    for i in 1:nbpart
        particles.case[i,1] = trunc(Int, particles.pos[i,1] / dimx * nx) + 1
        particles.case[i,2] = trunc(Int, particles.pos[i,2] / dimy * ny) + 1
    end


    # surface(compute_rho(particles, mesh))
    # savefig("rho.png")
    

    @show sum(compute_rho(particles, mesh)) 
    @show sum(particles.p) 

    @test maximum(abs.(rho .- compute_rho(particles, mesh))) < 1e-6
    
    
end
