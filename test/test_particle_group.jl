@testset "Particle Group 2D2V" begin

    n_particles = 10
    charge = -1.0
    mass = 1.0
    n_weights = 1

    particle_group =
        ParticleGroup{2,2}(n_particles, charge = 1.0, mass = 1.0, n_weights = 1)

    for ipart = 1:n_particles
        GEMPIC.set_x!(particle_group, ipart, [ipart, 0.0])
        GEMPIC.set_v!(particle_group, ipart, [ipart^2, 0.0])
        GEMPIC.set_weights!(particle_group, ipart, ipart / n_particles)
    end

    alpha = 0.1
    kx = 0.5

    mesh = TwoDGrid(0, 2π / kx, 128, 0, 1, 16)

    sampler = LandauDamping(alpha, kx)

    ParticleInCell.sample!(particle_group, mesh, sampler)

    @test true

end
