@testset "Poisson 2D" begin


eta1_min = .0; eta1_max = 2π
eta2_min = .0; eta2_max = 2π

nc_eta1 = 128; nc_eta2 = 128

ex = zeros(nc_eta1+1,nc_eta2+1)
ey = zeros(nc_eta1+1,nc_eta2+1)
ex_exact = zeros(nc_eta1+1,nc_eta2+1)
ey_exact = zeros(nc_eta1+1,nc_eta2+1)
rhs = zeros(nc_eta1+1,nc_eta2+1)
rho = zeros(nc_eta1+1,nc_eta2+1)
phi = zeros(nc_eta1+1,nc_eta2+1)
phi_exact = zeros(nc_eta1+1,nc_eta2+1)

println(" eta1_min, eta1_max, nc_eta1 ", eta1_min, eta1_max, nc_eta1)
println(" eta2_min, eta2_max, nc_eta2 ", eta2_min, eta2_max, nc_eta2)

grid = TwoDGrid( eta1_min, eta1_max, nc_eta1, 
                 eta2_min, eta2_max, nc_eta2) 

poisson = TwoDPoissonPeriodic( grid )

mode = 2
for i = 1:nc_eta1+1
   for j = 1:nc_eta2+1
      x1 = (i-1)*(eta1_max-eta1_min)/nc_eta1
      x2 = (j-1)*(eta2_max-eta2_min)/nc_eta2
      phi_exact[i,j] = real(mode,f64) * sin(mode*x1) * cos(mode*x2)
      ex_exact[i,j]  =  1._f64*real(mode,f64)**2*cos(mode*x1)*cos(mode*x2)
      ey_exact[i,j]  = -1._f64*real(mode,f64)**2*sin(mode*x1)*sin(mode*x2)
      rho[i,j] = -2._f64 * real(mode,f64)**3 * sin(mode*x1)*cos(mode*x2)
   end
end

rhs = rho
solve( poisson, phi, rhs)
error =  maximum(abs(phi_exact+phi))
rhs = rho
solve( poisson, phi, rhs)
error = maxval(abs(phi_exact+phi))
rhs = rho
solve( poisson, ex, ey, rhs)
error = maxval(abs(ex_exact-ex))
error = maxval(abs(ey_exact-ey))
rhs = rho
solve( poisson, ex, ey, rhs)
error = maxval(abs(ex_exact-ex))
error = maxval(abs(ey_exact-ey))

x1_min = 0._f64
x1_max = 1._f64

x2_min = 0._f64
x2_max = 1._f64

Nc_x1 = 32
Nc_x2 = 64

phi = zeros(Nc_x1+1,Nc_x2+1)
E1 = zeros(Nc_x1+1,Nc_x2+1)
E2 = zeros(Nc_x1+1,Nc_x2+1)
rho = zeros(Nc_x1+1,Nc_x2+1)

rho = 1

poisson = Poisson2DPeriodic(grid)

compute_phi_from_rho( phi, rho )

compute_e_from_rho( E1, E2, rho )

  nx, ny = 32, 32
  xmin, xmax = 0.0, 2π
  ymin, ymax = 0.0, 2π
  degree_smoother = 3
  n_particles = 1000

  grid = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)
  poisson = TwoDPoissonPeriodic( grid )
  
  # Initialize the kernel smoother
  kernel_smoother = ParticleMeshCoupling2D( grid, n_particles, degree_smoother, :collocation)

  @test true
  
  # Initialize the PIC field solver
  pic_poisson = ( poisson, num_cells, kernel_smoother)

  @test true

end
