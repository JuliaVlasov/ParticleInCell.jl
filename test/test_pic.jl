using ParticleInCell


nstep = 250
time = zeros(nstep + 1)
energy = zeros(nstep + 1)
dimx, dimy = 2π, 1.0
nx, ny = 128, 16
nbpart = 100 * nx * ny
dx, dy = dimx / nx, dimy / ny

# particles = landau( dimx, nbpart)

@time vm2d2v(1, time, energy)

nbpart = 10
particles = zeros(7, nbpart)
fields = ones(3, nx, ny)
interpolation!(particles, fields, nbpart, nx, ny, dx, dy)
