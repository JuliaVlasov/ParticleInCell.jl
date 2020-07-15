using Plots
using ParticleInCell

@testset "Particles" begin

    nx   = 10	
    ny   = 10  
    
    dimx = 10
    dimy = 10  
    
    mesh = Mesh( dimx, nx, dimy, ny)
    
    nbpart = 10
    particles = Particles(nbpart)

    particles.pos[:,1] .= LinRange(0,dimx,nbpart+1)[1:end-1] .+ 0.1
    particles.pos[:,2] .= LinRange(0,dimy,nbpart+1)[1:end-1] .+ 0.1

    update_cells!(particles, mesh)
    
    for i in 1:particles.nbpart

        @test particles.case[i,1] == i
        @test particles.case[i,2] == i
    
    end
    
end
