# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.2
#     language: julia
#     name: julia-1.10
# ---

# # Landau damping

# +
using ParticleInCell

function run()

    nx = 128
    ny = 16
    alpha = 0.1
    kx = 0.5
    ky = 0.0
    dimx = 2pi / kx
    dimy = 1.0

    mesh = TwoDGrid(dimx, nx, dimy, ny)
   
    ex = zeros(nx + 1, ny + 1)
    ey = zeros(nx + 1, ny + 1)
    bz = zeros(nx + 1, ny + 1)
    jx = zeros(nx + 1, ny + 1)
    jy = zeros(nx + 1, ny + 1)

    nbpart = 100 * nx * ny

    group = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    alpha = 0.1
    kx = 0.5

    sampler = LandauDamping(alpha, kx)

    ParticleInCell.sample!(group, mesh, sampler)

    fdtd = FDTD(mesh)

    time = 0
    dt = 0.01

    for i = 1:nx, j = 1:ny+1
        fdtd.ex[i, j] = alpha / kx * sin(kx * mesh.x[i])
    end

    kernel = CloudInCell()

    particles = group.array

    for istep = 1:100

        istep > 1 && faraday!(fdtd, mesh, 0.5dt)

        update_fields!(ex, ey, bz, mesh, fdtd)

        push_v!(group, kernel, mesh, ex, ey, bz, dt)

        push_x!(group, mesh, 0.5dt)

        compute_current!(jx, jy, mesh, kernel, group)

        push_x!(group, mesh, 0.5dt)

        faraday!(fdtd, mesh, 0.5dt)

        ampere_maxwell!(fdtd, mesh, jx, jy, dt)

        time = time + dt

    end

end

# -

run()


