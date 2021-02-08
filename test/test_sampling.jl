using Pkg
Pkg.activate("/Users/navaro/JuliaProjects/ParticleInCell.jl")

using Plots

using Revise

using ParticleInCell

nx        = 256	    # nombre de pts suivant x
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

ex  = zeros(nx,ny)
rho = zeros(nx,ny)

for i=1:nx
    aux1 = alpha/kx * sin(kx*mesh.x[i])
    aux2 = alpha * cos(kx*mesh.x[i])
    for j=1:ny
        ex[i,j] = aux1
        rho[i,j] = aux2
    end
end

@show sum(rho[1:nx,1:ny])

nbpart = 500*nx*ny

particles = Particles(nbpart)

landau_sampling!( particles, alpha, kx )

histogram(particles.vit[:,1], normalize=true)

update_cells!( particles, mesh)


surface(compute_rho(particles, mesh))

surface(rho)


