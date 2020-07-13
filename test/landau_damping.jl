using OffsetArrays
using Sobol
using Plots
import ParticleInCell: Mesh, Fields, Particles, landau_sampling!
import ParticleInCell: decalage!, faraday!, ampere!, interpol_eb!
import ParticleInCell: push_x!, push_v!

@testset "Landau Damping" begin

    nx        = 128	    # nombre de pts suivant x
    ny        = 16   	# nombre de pts suivant y
    cfl       = 0.9 	# nombre de Courant-Friedrich-Levy
    tfinal    = 50     	# temps final
    nstepmax  = 2000000	# nbre d'iterations maxi
    exext     = 0		# champ electrique exterieur suivant x
    eyext     = 0	 	# champ electrique exterieur suivant y
    bzext     = 0		# champ magnetique exterieur
    charge    = 1       # charge d'une macro particule
    masse     = 1       # masse d'une macro particule
    c         = 8	    # vitesse de la lumiere
    e0        = 1       # permittivite du vide
    
    csq = c*c
    q_over_m = charge / masse
    poids = charge
    
    alpha = 0.1
    kx = 0.5
    ky = 0.
    dimx = 2*pi/kx
    dimy = 1  
    poids = dimx * dimy 
    
    mesh = Mesh( dimx, nx, dimy, ny)
    
    dx = maximum(mesh.hx)
    dy = maximum(mesh.hy)
    dt = cfl  / sqrt(1/(dx*dx)+1/(dy*dy)) / c
    
    nstep = trunc(Int, tfinal ÷ dt)
    
    println(" cfl = $cfl ")
    println(" dx = $dx dy = $dy dt = $dt")
    println(" Nnumber of steps $nstep ")
    
    fields = Fields(nx, ny)
    
    ex = OffsetArray{Float64}(undef, 0:nx,0:ny)
    ey = OffsetArray{Float64}(undef, 0:nx,0:ny)
    bz = OffsetArray{Float64}(undef, 0:nx,0:ny)
    jx = OffsetArray{Float64}(undef, 0:nx,0:ny)
    jy = OffsetArray{Float64}(undef, 0:nx,0:ny)
    
    time  = 0
    iplot = 0
    
    nstep > nstepmax && ( nstep = nstepmax )
    
    istep = 1
    
    fill!(fields.ex, 0)
    fill!(fields.ey, 0)
    fill!(fields.bz, 0)
    fill!(fields.jx, 0)
    fill!(fields.jy, 0)
    fill!(fields.r0, 0)
    fill!(fields.r1, 0)
    
    for i=0:nx-1
        aux1 = alpha/kx * sin(kx*mesh.x[i])
        aux2 = alpha * cos(kx*mesh.x[i])
        for j=0:ny
            fields.ex[i,j] = aux1
            fields.r1[i,j] = aux2
        end
    end
          
    nbpart = 100*nx*ny
    
    particles = Particles(nbpart)
    
    landau_sampling!( particles, alpha, kx )
    
    xp = view(particles.pos, :, 1)
    vp = particles.vit
    
    pp = plot(layout=(3,1))
    histogram!(pp[1,1], xp, normalize=true, bins = 100, lab="x")
    plot!(pp[1,1], x -> (1+alpha*cos(kx*x))/(2π/kx), 0., 2π/kx, lab="")
    histogram!(pp[2,1], vp[:,1], normalize=true, bins = 100, lab="vx")
    plot!(pp[2,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
    histogram!(pp[3,1], vp[:,2], normalize=true, bins = 100, lab="vy")
    plot!(pp[3,1], v -> exp( - v^2 / 2) * 4 / π^2 , -6, 6, lab="")
    savefig("particles.png")
    
    
    for istep = 1:nstep
    
       if istep > 1
           faraday!( fields, 0.5dt ) 	#Calcul de B(n-1/2) --> B(n)			
       end
    
       decalage!( mesh, fields, ex, ey, bz )
       interpol_eb!( ex, ey, bz, particles, mesh )
    
       push_v!( particles )
    
       psu_x!( particles, 0.5)  # x(n) --> x(n+1/2)
       output_particles!( particles )
       calcul_j_cic( particles, fields, jx, jy )
       push_x!( particles, 0.5)  # x(n+1/2) -- x(n+1)
       output_particles!( particles )
            
       faraday!( mesh, fields, dt, c, e0 )   #Calcul de B(n) --> B(n+1/2)
       ampere!( mesh, fields, dt, c, e0 )    #Calcul de E(n) --> E(n+1)
       boundary_conditions!( mesh, fields )
    
       time = time + dt
    
    
    end 

end
