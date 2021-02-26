using Test
using ParticleInCell

@testset "ParticleMeshCoupling 2D" begin


  n_cells = 10 # Number of cells
  n_particles = 4 # Number of particles
  spline_degree = 3 # Spline degree
  xmin, xmax, nx = 0.0, 2.0, n_cells
  ymin, ymax, ny = 0.0, 1.0, n_cells

  mesh = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)

  volume = (xmax - xmin) * (ymax - ymin)

  x_vec = [0.1 0.65 0.7 1.5; 0.0 0.0 0.0 0.0]
  v_vec = [1.5 0.00 0.0 0.0; 0.0 0.5 0.0 0.0]

  pg = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)

  @show size(pg.array)
  
  for i_part = 1:n_particles
     pg.array[1, i_part] = x_vec[1,i_part]
     pg.array[2, i_part] = x_vec[2,i_part]
     pg.array[3, i_part] = v_vec[1,i_part]
     pg.array[4, i_part] = v_vec[2,i_part]
     pg.array[5, i_part] = 1.0 / n_particles
  end 

  # Initialize the kernel
  kernel = ParticleMeshCoupling2D( pg, mesh, spline_degree, :collocation)

  # Reference values of the shape factors
  index_grid = zeros(Int, (2, 4))
  index_grid[1,:] .= [-2, 1, 1, 5]
  index_grid[2,:] .= [-3, -3, -3, -3]
  

  values_grid = zeros(Float64,(4,2,4))
  values_grid[:,1,1] .= [ 2.0833333333333332E-002, 0.47916666666666663, 
                          0.47916666666666663, 2.0833333333333332E-002]
  values_grid[:,1,3] .= values_grid[:,1,1]
  values_grid[:,1,4] .= values_grid[:,1,1] 
  values_grid[:,1,2] .= [ 7.0312500000000000E-002, 0.61197916666666663, 
                         0.31510416666666663, 2.6041666666666665E-003]
  values_grid[1,2,:] .= 0.0
  values_grid[2,2,:] .= 1.0/6.0
  values_grid[3,2,:] .= 2.0/3.0
  values_grid[4,2,:] .= 1.0/6.0
  
  @show values_grid

  rho_dofs = zeros(Float64, 100)
  rho_dofs1 = zeros(Float64, 100)
  rho_dofs_ref = zeros(Float64, 100)
  rho_dofs_pp = zeros(Float64, (16,100))

  # Accumulate rho
  for i_part in 1:n_particles
     x = pg.array[1,i_part]
     y = pg.array[2,i_part]
     w = pg.array[5,i_part]
     add_charge!(rho_dofs,  kernel, x, y, w)
     # add_charge_pp!(rho_dofs1, kernel, x, y, w)
  end

  rho_dofs_ref[8:10] .= values_grid[1:3,1,1]
  rho_dofs_ref[1]     = values_grid[4,1,1]
  rho_dofs_ref[1:4]  .= rho_dofs_ref[1:4] .+ values_grid[:,1,2] .+ values_grid[:,1,3]
  rho_dofs_ref[5:8]  .= rho_dofs_ref[5:8] .+ values_grid[:,1,4]

  rho_dofs_ref[71:80]  .= rho_dofs_ref[1:10] ./ 6.0
  rho_dofs_ref[81:90]  .= rho_dofs_ref[1:10] .* 2.0/3.0
  rho_dofs_ref[91:100] .= rho_dofs_ref[1:10] ./ 6.0
  rho_dofs_ref[1:10]   .= 0.0

  rho_dofs_ref .*= n_cells^2 / volume / n_particles

  @test maximum(abs.(rho_dofs  .- rho_dofs_ref)) < 1e-14
  # @test maximum(abs.(rho_dofs1 .- rho_dofs_ref)) < 1e-14 

#=
  spline_pp_b_to_pp_2d(kernel.spline_pp,[n_cells,n_cells],rho_dofs,rho_dofs_pp)
  
  # Test function evaluation
  for i_part = 1:n_particles
     xi = get_x(particle_group, i_part)
     evaluate( kernel, xi, rho_dofs, particle_values[i_part])
     evaluate_pp(kernel, xi, rho_dofs_pp, particle_values1[i_part])
  end
  particle_values_ref = [1.1560058593749998, 2.3149278428819446,
                         2.2656250000000000, 1.1512586805555554] ./ volume

  particle_values_ref = 0.0
  for i_part = 1:n_particles
     for i=1:4
        i1 = mod(index_grid[1,i_part]+i-2, n_cells)
        for j=1:4
           i2 = mod(index_grid[2,i_part]+j-2, n_cells)
           particle_values_ref[i_part] += (values_grid[i,1,i_part]
                                          *values_grid[j,2,i_part]
                                          *rho_dofs_ref[i1+i2*n_cells+1])
        end
     end
  end

  error  = maximum(abs.(particle_values  .- particle_values_ref))
  error1 = maximum(abs.(particle_values1 .- particle_values_ref))

  @test error  < 1.e-14
  @test error1 < 1.e-14

=#
end 
