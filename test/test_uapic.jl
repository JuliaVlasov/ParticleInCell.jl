using UAPIC
using Test

@testset " Test particles input file read " begin


mesh = Mesh( 0, 4π, 128, 0, 2π, 64 )

particles = read_particles("particles.dat", mesh)

@test particles.nbpart == 204800


end

@testset " Test Particles-MeshFields interaction " begin

     xmin, xmax = 0.0, 20.0
     ymin, ymax = 0.0, 20.0
     nx,   ny   = 20, 20

     mesh = Mesh( xmin, xmax, nx, ymin, ymax, ny )

     dx, dy = mesh.dx, mesh.dy

     fields = MeshFields( mesh )

     nbpart = 121
     w      = 1/nbpart

     particles = Particles( nbpart, w )

     k = 1
     for i = 5:nx-5, j = 5:ny-5

         particles.x[1,k] = (i-0.5)*dx
         particles.x[2,k] = (j-0.5)*dx

         k += 1

     end

     compute_rho_m6!( fields, particles )

     @test integrate( fields.ρ, mesh ) ≈ 0.0 atol = 1e-4

     for i = 1:nx+1, j = 1:ny+1

         fields.e[1,i,j] = (i-1)*dx
         fields.e[2,i,j] = (j-1)*dy

     end
     
     interpol_eb_m6!( particles, fields )

     err_x, err_y = 0.0, 0.0

     for k = 1:nbpart
         
         xp = particles.x[1,k]
         yp = particles.x[2,k]
         ex = particles.e[1,k]
         ey = particles.e[2,k]

         err_x += abs( particles.e[1,k] - xp ) 
         err_y += abs( particles.e[2,k] - yp ) 

     end

     err_x /= nbpart
     err_y /= nbpart

     @test err_x ≈ 0.0 atol = 1e-6
     @test err_y ≈ 0.0 atol = 1e-6

end
