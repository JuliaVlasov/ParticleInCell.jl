@testset "Particle Group 2D2V" begin

    n_particles = 10
    charge = -1.0
    mass = 1.0
    n_weights = 1

    particle_group =
        ParticleGroup{2,2}(n_particles, charge = 1.0, mass = 1.0, n_weights = 1)

    for ipart = 1:n_particles
        ParticleInCell.set_x!(particle_group, ipart, [ipart, 0.0])
        ParticleInCell.set_v!(particle_group, ipart, [ipart^2, 0.0])
        ParticleInCell.set_w!(particle_group, ipart, ipart / n_particles)
    end

    alpha = 0.1
    kx = 0.5

    sampler = LandauDamping(alpha, kx)
    sample!(particle_group, sampler)

    @test true

end
