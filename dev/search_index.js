var documenterSearchIndex = {"docs":
[{"location":"vlasov-maxwell/#Vlasov-Maxwell","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell","text":"","category":"section"},{"location":"vlasov-maxwell/#Landau-damping","page":"Vlasov-Maxwell 2D","title":"Landau damping","text":"","category":"section"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"using Plots, LinearAlgebra","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"using ParticleInCell\n\nnx = 128 \nny = 16   \n\nalpha = 0.1\nkx = 0.5\nky = 0.\ndimx = 2*pi/kx\ndimy = 1  \npoids = dimx * dimy \n\nmesh = TwoDGrid( dimx, nx, dimy, ny)\nfdtd = FDTD(mesh)\n\ntime  = 0\n\nfor i=1:nx, j=1:ny\n    fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)\nend\nsurface(fdtd.ex )","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"nbpart = 100*nx*ny\nparticles = zeros(7, nbpart)\nlandau_sampling!( particles, alpha, kx )\n\np = plot(layout=4)\nhistogram!(p[1], particles[1,:], normalize=true, label=\"x\")\nhistogram!(p[2], particles[2,:], normalize=true, label=\"y\")\nhistogram!(p[3], particles[3,:], normalize=true, label=\"vx\")\nhistogram!(p[4], particles[4,:], normalize=true, label=\"vy\")","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"compute_current!( mesh, particles)\n\np = plot(layout=2)\nsurface!(p[1], mesh.jx)\nsurface!(p[2], mesh.jy)","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"function run( fdtd, particles, mesh, nstep, dt)\n    \n    alpha = 0.1\n    kx = 0.5\n    landau_sampling!( particles, alpha, kx )\n    for i=1:nx, j=1:ny\n        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)\n        fdtd.ey[i,j] = 0.0\n        fdtd.bz[i,j] = 0.0\n    end\n\n    time = 0\n    energy = Float64[0.5 * log( sum( fdtd.ex.^2) * mesh.dx * mesh.dy)]\n    t = Float64[time]\n    \n    for istep in 1:nstep\n    \n       istep > 1 && faraday!( fdtd, mesh, 0.5dt ) \n       interpolation!( particles, mesh )\n       push_v!( particles, dt )\n       push_x!( particles, mesh, 0.5dt) \n       compute_current!( mesh, particles)\n       push_x!( particles, mesh, 0.5dt) \n       faraday!(fdtd, mesh, 0.5dt)\n       ampere_maxwell!(fdtd, mesh, dt)\n       time = time + dt\n       push!(t, time)\n       push!(energy, 0.5 * log( sum(fdtd.ex.^2) * mesh.dx * mesh.dy))\n    \n    end\n   \n    t, energy\n    \nend","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"dt = 0.01\nnstep = 1000\nt, energy = run( fdtd, particles, mesh, nstep, dt)\nplot(t, energy)","category":"page"},{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [ParticleInCell]\nOrder   = [:function]","category":"page"},{"location":"functions/#ParticleInCell.add_charge!-Tuple{Any,ParticleMeshCoupling2D,Any,Any,Any}","page":"Functions","title":"ParticleInCell.add_charge!","text":"add_charge(pm, x_position, y_position wp, ρ_dofs)\n\nAdd charge of single particle\n\nposition : Particle position\nwp : Particle weight times charge\nρ_dofs : spline coefficient of accumulated density\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_charge!-Tuple{Array{Float64,1},ParticleMeshCoupling1D,Float64,Float64}","page":"Functions","title":"ParticleInCell.add_charge!","text":"add_charge!( rho, p, position, marker_charge)\n\nAdd charge of one particle\n\np             : kernel smoother object\nposition      : Position of the particle\nmarker_charge : Particle weights time charge\nrho_dofs      : Coefficient vector of the charge distribution\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_charge!-Tuple{PICPoisson2D,Any,Any}","page":"Functions","title":"ParticleInCell.add_charge!","text":"add_charge!(pic, position, marker_charge)\n\nAdd charge from one particle\n\nself : Pic Poisson solver object\nposition : Position of the particle\nmarker_charge : Particle weight times charge\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_charge_pp!-Tuple{Any,ParticleMeshCoupling2D,Any,Any,Any}","page":"Functions","title":"ParticleInCell.add_charge_pp!","text":"add_charge_single_spline_pp_2d(pm, position, wp, ρ_dofs)\n\nAdd charge of single particle\n\nInformation about the 2d mesh\ndelta_x(2)  : Value of grid spacing along both directions.\ndomain(2,2) : Definition of the domain: domain(1,1) = x1min, domain(2,1) = x2min,  domain(1,2) = x1max, domain(2,2) = x2max\nInformation about the particles\nno_particles : Number of particles of underlying PIC method (processor local)\nn_span : Number of intervals where spline non zero (degree + 1)\nscaling\nposition : Particle position\nwp : Particle weight times charge\nρ_dofs : spline coefficient of accumulated density\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_charge_pp!-Tuple{Array{Float64,1},ParticleMeshCoupling1D,Float64,Any}","page":"Functions","title":"ParticleInCell.add_charge_pp!","text":"add_charge_pp!(rho_dofs, p, position, marker_charge)\n\nAdd charge of one particle\n\np             : kernel smoother object\nposition      : Position of the particle\nmarker_charge : Particle weights time charge\nrho_dofs      : Coefficient vector of the charge distribution\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_current_update_v!-Tuple{AbstractArray,ParticleMeshCoupling1D,Float64,Float64,Float64,Float64,Array{Float64,1},Float64}","page":"Functions","title":"ParticleInCell.add_current_update_v!","text":"add_current_update_v!( j_dofs, p, \n                       position_old, position_new, \n                       marker_charge, qoverm, \n                       bfield_dofs, vi)\n\nAdd current for one particle and update v (according to H_p1 part in Hamiltonian splitting)\n\nRead out particle position and velocity\nCompute index_old, the index of the last DoF on the grid the \n\nparticle contributes to, and r_old, its position (normalized to cell size one).\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.add_current_update_v_pp!-Tuple{AbstractArray,ParticleMeshCoupling1D,Any,Any,Float64,Float64,Array{Float64,1},Array{Float64,1}}","page":"Functions","title":"ParticleInCell.add_current_update_v_pp!","text":"add_current_update_v_pp!( j_dofs, p, position_old, position_new, \n                          marker_charge, qoverm, bfield_dofs, vi)\n\nAdd current for one particle and update v  (according to H_p1 part in Hamiltonian splitting)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.b_to_pp-Tuple{ParticleInCell.SplinePP,Int64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.b_to_pp","text":"b_to_pp( SplinePP, ncells, b_coeffs)\n\nConvert 1d spline in B form to spline in pp form with  periodic boundary conditions\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.b_to_pp_1d_cell!-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.b_to_pp_1d_cell!","text":"b_to_pp_1d_cell( self, b_coeffs, pp_coeffs )\n\nConvert 1d spline in B form in a cell to spline in pp form with periodic boundary conditions\n\nspline : arbitrary degree 1d spline \nb_coeffs(self%degree+1) : coefficients of spline in B-form\npp_coeffs(self%degree+1) : coefficients of spline in pp-form\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.b_to_pp_2d!-Tuple{Any,ParticleInCell.SplinePP,ParticleInCell.SplinePP,Any}","page":"Functions","title":"ParticleInCell.b_to_pp_2d!","text":"b_to_pp_2d!( pp, spl1, spl2, b)\n\nConvert 2d spline in B form to spline in pp form   \n\nn_cells(2) : number of gridcells\nb_coeffs   : coefficients of spline in B-form\npp_coeffs  : coefficients of spline in pp-form\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.b_to_pp_2d_cell!-Tuple{Any,ParticleInCell.SplinePP,ParticleInCell.SplinePP,Any,Any,Any}","page":"Functions","title":"ParticleInCell.b_to_pp_2d_cell!","text":"b_to_pp_2d_cell(spline1, spline2, b_coeffs, pp_coeffs, i, j)\n\nConvert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions \n\nspline1 : arbitrary degree 1d spline\nspline2 : arbitrary degree 1d spline\nn_cells(2) : number of gridcells\nbcoeffs(ncells(1)*n_cells(2)) : coefficients of spline in B-form\nppcoeffs((spline1.degree+1)*(spline2.degree+1),ncells(1)*n_cells(2)) : coefficients of spline in pp-form\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_b_from_e!-Tuple{Array{Float64,1},Maxwell1DFEM,Float64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.compute_b_from_e!","text":"compute_b_from_e!( field_out, maxwell_solver, delta_t, field_in)\n\nCompute Bz from Ey using strong 1D Faraday equation for spline coefficients\n\nB_z^new(x_j) = B_z^old(x_j) - fracDelta tDelta x (E_y(x_j) - E_y(x_j-1)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_e_from_b!-Tuple{Array{Float64,1},Maxwell1DFEM,Float64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.compute_e_from_b!","text":"compute_e_from_b!(field_out, maxwell_solver, delta_t, field_in)\n\ncompute Ey from Bz using weak Ampere formulation \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_e_from_j!-Tuple{Array{Float64,1},Maxwell1DFEM,Array{Float64,1},Int64}","page":"Functions","title":"ParticleInCell.compute_e_from_j!","text":"compute_e_from_j!(e, maxwell_solver, current, component)\n\nCompute E_i from j_i integrated over the time interval using weak Ampere formulation\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_rhs_from_function!-Tuple{Array{Float64,1},Maxwell1DFEM,Function,Int64}","page":"Functions","title":"ParticleInCell.compute_rhs_from_function!","text":"computerhsfromfunction(self, func, degree, coefsdofs)\n\nCompute the FEM right-hand-side for a given function f and periodic splines of given degree.\n\nIts components are int f N_i dx where N_i is the B-spline starting at x_i. \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_shape_factor-Tuple{ParticleMeshCoupling2D,Any,Any}","page":"Functions","title":"ParticleInCell.compute_shape_factor","text":"computeshapefactor(pm, xp, yp)\n\nHelper function computing shape factor\n\npm : kernel smoother object\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.eval_uniform_periodic_spline_curve-Tuple{Int64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.eval_uniform_periodic_spline_curve","text":"eval_uniform_periodic_spline_curve( degree, scoef )\n\nEvaluate uniform periodic spline curve defined by coefficients scoef at  knots (which are the grid points) \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate-NTuple{4,Any}","page":"Functions","title":"ParticleInCell.evaluate","text":"evaluate_field_single_spline_2d(pm, position, field_dofs)\n\nposition(pm.dim) : Position where to evaluate\nfielddofs(pm.ndofs) : Degrees of freedom in kernel representation.\nfield_value : Value of the field\n\nEvaluate field with given dofs at position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate-Tuple{ParticleMeshCoupling1D,Float64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.evaluate","text":"evaluate(p, position, field_dofs)\n\nEvaluate field at position\n\np : Kernel smoother object \nposition : Position of the particle\nfield_dofs : Coefficient vector for the field DoFs\nfield_value : Value(s) of the electric fields at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_fields-Tuple{Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_fields","text":"evaluate_field_single_2d(pic, position, components, func_value)\n\nEvaluate components of the electric field as one position\n\nself : PIC Poisson solver object\ncomponents\nposition : Position of the particle\nfunc_value\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_multiple-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_multiple","text":"evaluate_multiple(pm, position, components, field_dofs)\n\nEvaluate multiple fields at position \u0007 position\n\nposition(pm%dim) : Position where to evaluate\ncomponents(:) : Components of the field that shall be evaluated\nfield_dofs(:,:) : Degrees of freedom in kernel representation.\nfield_value(:) : Value of the field\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_pp-Tuple{ParticleMeshCoupling1D,Float64,Array{Float64,2}}","page":"Functions","title":"ParticleInCell.evaluate_pp","text":"evaluate_pp(p, position, field_dofs_pp)\n\nEvaluate field at position using horner scheme\n\np : Kernel smoother object \nposition : Position of the particle\nfield_dofs_pp : Degrees of freedom in kernel representation.\nfield_value : Value(s) of the electric fields at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_pp-Tuple{ParticleMeshCoupling2D,Any,Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_pp","text":"evaluate_pp(pm, position, field_dofs_pp)\n\nEvaluate field at position using horner scheme\n\npm : kernel smoother object    \nposition : Position where to evaluate\nfielddofspp(:,:) : Degrees of freedom in kernel representation.\nfield_value : Value of the field\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_rho!-Tuple{Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_rho!","text":"evaluate_rho!(pic, position)\n\nEvaluate charge density at rho at one position\n\nself : Pic Poisson solver object\nposition : Position of the particle\nfunc_value : Value of rho at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_rho-Tuple{Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_rho","text":"evaluate_phi!(pic, position)\n\nEvaluate potential at one position\n\nself : Pic Poisson solver object\nposition : Position of the particle\nfunc_value : Value of phi at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_cell_and_offset-Tuple{OneDGrid,Any}","page":"Functions","title":"ParticleInCell.get_cell_and_offset","text":"get_cell_and_offset( mesh, x )\n\nGet cell and offset\n\nWe compute the cell indices where the particle is and its relative  normalized position inside the cell\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_charge-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_charge","text":"get_charge( p, i; i_wi=1)\n\nGet charge of ith particle of p (q * particle_weight)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_mass-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_mass","text":"get_mass( p, i; i_wi=1)\n\nGet mass of ith particle of p (m * particle_weight)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_spin-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_spin","text":"get_spin( p, i, j)\n\nGet the jth weight of the ith particle weights of group p\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_v-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_v","text":"get_v( p, i )\n\nGet velocity of ith particle of p\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_weights-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_weights","text":"get_weights( p, i)\n\nGet ith particle weights of group p\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_x-Tuple{OneDGrid,Any}","page":"Functions","title":"ParticleInCell.get_x","text":"get_x( mesh, i )\n\nGet position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.get_x-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64}} where V where D","page":"Functions","title":"ParticleInCell.get_x","text":"get_x( p, i )\n\nGet position of ith particle of p\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.horner_1d-Tuple{Int64,Any,Float64,Int64}","page":"Functions","title":"ParticleInCell.horner_1d","text":"horner_1d(degree, pp_coeffs, x, index)\n\nPerform a 1d Horner schema on the pp_coeffs at index\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.horner_2d-NTuple{5,Any}","page":"Functions","title":"ParticleInCell.horner_2d","text":"horner_2d(degrees, pp_coeffs, position, indices, ncells)\n\nPerform a 2d hornerschema on the pp_coeffs at the indices\n\ndegree : degree of the spline\npp_coeffs : coefficients of spline in pp-form\nposition(2) : point at which we evaluate our spline\nindices(2) : indices of cell in which is x\nncells(2) : number of gridcells\nres : value of the splinefunction at position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.horner_m_2d!-Tuple{Any,ParticleInCell.SplinePP,ParticleInCell.SplinePP,Any,Any,Any}","page":"Functions","title":"ParticleInCell.horner_m_2d!","text":"horner_m_2d!(val, spl1, spl2, degree, x)\n\nPerform two times a 1d hornerschema on the poly_coeffs\n\nval : array of values\ndegree : degree of the spline\nx : point at which we evaluate our spline\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.horner_primitive_1d-Tuple{Array{Float64,1},Any,Any,Any}","page":"Functions","title":"ParticleInCell.horner_primitive_1d","text":"horner_primitive_1d(val, degree, pp_coeffs, x)\n\nPerform a 1d Horner schema on the pp_coeffs evaluate at x\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.index_1dto2d_column_major-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.index_1dto2d_column_major","text":"index_1dto2d_column_major(pm, index1d)\n\nSelf function computes the index of a 1D array that stores 2D data in column major ordering.  It also takes periodic boundary conditions into account.\n\nindex1d_1 !< indice along x (start counting with zero).\nindex1d_2 !< indice along y (start counting with zero).\nindex2d   !< Corresponding index in 1d array representing 2d data (start counting with one).\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.inner_product-NTuple{4,Any}","page":"Functions","title":"ParticleInCell.inner_product","text":"inner_product( maxwell_solver, coefs1_dofs, coefs2_dofs, degree )\n\nmaxwell_solver : Maxwell solver object\ncoefs1_dofs : Coefficient for each DoF\ncoefs2_dofs : Coefficient for each DoF\n`degree : Specify the degree of the basis functions\n\nreturn squared L2 norm\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.l2norm_squared-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.l2norm_squared","text":"l2norm_squared(maxwell_solver, coefs_dofs, degree)\n\nCompute square of the L2norm \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.l2norm_squared2-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.l2norm_squared2","text":"l2norm_squared(maxwell_solver, coefs_dofs, degree)\n\nCompute square of the L2norm \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.l2projection!-Tuple{Array{Float64,1},Maxwell1DFEM,Function,Int64}","page":"Functions","title":"ParticleInCell.l2projection!","text":"l2projection!(coefs_dofs, maxwell, func, degree)\n\nCompute the L2 projection of a given function f on periodic splines  of given degree\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.operator_t!-Tuple{SplittingOperator,Any}","page":"Functions","title":"ParticleInCell.operator_t!","text":"operator_t(split, dt)\n\nPush x \n\nsplit :: time splitting object \ndt   :: time step\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_spin-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Int64,Any}} where V where D","page":"Functions","title":"ParticleInCell.set_spin","text":"set_spin( p, i, j)\n\nSet the jth weight of the ith particle weights of group p\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_v-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Array{Float64,1}}} where V where D","page":"Functions","title":"ParticleInCell.set_v","text":"set_v( p, i, v)\n\nSet velocity of ith particle of p to v\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_v-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Float64}} where V where D","page":"Functions","title":"ParticleInCell.set_v","text":"set_v( p, i, v)\n\nSet velocity of ith particle of p to v\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_weights-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Array{Float64,1}}} where V where D","page":"Functions","title":"ParticleInCell.set_weights","text":"set_weights( p, i, w)\n\nSet weights of ith particle of p to w\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_weights-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Float64}} where V where D","page":"Functions","title":"ParticleInCell.set_weights","text":"set_weights( p, i, w)\n\nSet weights of particle @ i\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_x-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Array{Float64,1}}} where V where D","page":"Functions","title":"ParticleInCell.set_x","text":"set_x( p, i, x )\n\nSet position of ith particle of p to x \n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.set_x-Union{Tuple{V}, Tuple{D}, Tuple{ParticleGroup{D,V},Int64,Float64}} where V where D","page":"Functions","title":"ParticleInCell.set_x","text":"set_x( p, i, x)\n\nSet position of ith particle of p to x\n\nnote: Note\nif x is a scalar value, only the first x dimension will be set.\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.solve!","text":"Compute electric field from charge density\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any,Any,Poisson2DPeriodic,Any}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( ex, ey, poisson, rho )\n\nsolves Poisson equation to compute electric fields\n\nE(xy) = -nabla phi(xy) \n-Delta phi(xy) = rho(xy)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any,Poisson2DPeriodic,Any}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( poisson, phi, rho )\n\ncomputes phi from rho \n\n-Delta phi(xy) = rho(xy)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( pic )\n\nSolve for phi and fields\n\npoisson : Pic Poisson solver object\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_fields!-Tuple{Any}","page":"Functions","title":"ParticleInCell.solve_fields!","text":"solve_fields!( pic )\n\nSolve efields from rho\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_fields!-Tuple{SplittingOperator}","page":"Functions","title":"ParticleInCell.solve_fields!","text":"Solve Poisson's equation for the electric field\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_phi!-Tuple{Any}","page":"Functions","title":"ParticleInCell.solve_phi!","text":"solve_phi!( pic )\n\nSolve for potential\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.strang_splitting-Tuple{SplittingOperator,Any}","page":"Functions","title":"ParticleInCell.strang_splitting","text":"Strang splitting\n\nsplit :: time splitting object \ndt   :: time step\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.uniform_bsplines_eval_basis-Tuple{Int64,Float64}","page":"Functions","title":"ParticleInCell.uniform_bsplines_eval_basis","text":"uniform_bsplines_eval_basis( spline_degree, normalized_offset, bspl )\n\nUNIFORM B-SPLINE FUNCTIONS\n\nEvaluate all non vanishing uniform B-Splines in unit cell.\n\nReturns an array with the values of the b-splines of the  requested degree, evaluated at a given cell offset. The cell size is normalized between 0 and 1, thus the offset given must be a number between 0 and 1.\n\nOutput: \n\nbspl(1d+1)= B_d(-(d+1)2+d+x)B_d(-(d+1)2+x)\n\nwith d=spline_degree and x=normalized_offset where B_d=B_d-1*B_0 and B_0=1_-1212 and * is convolution the following FORTRAN code can be used for comparison with  deboor\n\ndo i=-d,d+1\n    t(i+d+1)=real(i,8)\nend do\ncall bsplvb(t,d+1,1,normalized_offset,d+1,out)\n\nWe also have the property (from the symmetry of the B-spline)\n\nout1d+1= B_d(-(d+1)2+xx)B_d(-(d+1)2+d+xx) \n\nwhere xx=1- normalized_offset\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.update_jv!-Tuple{AbstractArray,ParticleMeshCoupling1D,Float64,Float64,Int64,Float64,Float64,Float64,Float64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.update_jv!","text":"update_jv!(j_dofs, p, \n           lower, upper, index, marker_charge, \n           qoverm, sign, vi, bfield_dofs)\n\nHelper function for add_current_update_v.\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.update_jv_pp!-Tuple{AbstractArray,ParticleMeshCoupling1D,Float64,Float64,Int64,Float64,Float64,Float64,Array{Float64,1}}","page":"Functions","title":"ParticleInCell.update_jv_pp!","text":"update_jv_pp!( j_dofs, p, lower, upper, index, marker_charge, \n               qoverm, vi, bfield_dofs)\n\nHelper function for add_current_update_v.\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.update_positions!-Tuple{Any,Any,Any}","page":"Functions","title":"ParticleInCell.update_positions!","text":"update particle position xp\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.update_velocities!-NTuple{4,Any}","page":"Functions","title":"ParticleInCell.update_velocities!","text":"update particle velocities vp\n\n\n\n\n\n","category":"method"},{"location":"tsi/#Two-stream-instability","page":"Two-stream instability","title":"Two-stream instability","text":"","category":"section"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"using Plots\nusing Random","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"using ParticleInCell\n\nconst dt = 0.005     # Time step\nconst nt = 10000     # Number of time steps\nconst L  = 20π       #  Domain size \nconst nx = 320       # Number of grid cells\nconst np = nx * 20   # Number of particles\n\n\nmesh = Mesh( 0, 20π, nx)\nrng = MersenneTwister(42)\npoisson = Poisson1D( mesh )\nparticles = tsi(rng, mesh, np )\npm = ParticleMeshCoupling(particles, mesh)","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"function main()\n\n    mesh = Mesh( 0, 20π, nx)\n    poisson = Poisson1D( mesh )\n    rng = MersenneTwister(42)\n    pa = tsi(rng, mesh, np )\n    pm = ParticleMeshCoupling(pa, mesh)\n    energy = Float64[]\n    e = zeros(Float64, nx)\n    ρ = zeros(Float64, nx)\n    xmin = mesh.xmin\n    xmax = mesh.xmax\n    \n    for it in 1:nt+1\n        \n        update_positions!(pa, mesh, dt)\n        mat = compute_coeffs(pm, pa)\n        compute_rho!(ρ, mat, mesh, pa)\n        solve!(e, poisson, ρ)\n        update_velocities!(pa, e, mat, dt)\n        push!(energy, 0.5 * sum(e.^2) * mesh.dx) \n\n    end\n\n    energy\n\nend","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"results = main()\nt = (0:nt) .* dt\nplot( t, results, yaxis=:log)","category":"page"},{"location":"contents/#Contents","page":"Contents","title":"Contents","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"contents/#Index","page":"Contents","title":"Index","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"landau_damping/#Landau-damping","page":"Landau damping","title":"Landau damping","text":"","category":"section"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"using Plots\nusing Random","category":"page"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"\nusing ParticleInCell\n\nfunction main(nt, dt)\n    \n    nx = 50\n    np = 10000 * nx\n    mesh = Mesh( 0, 4π, nx)\n    poisson = Poisson1D( mesh )\n    rng = MersenneTwister(42)\n    α = 0.5\n    kx = 0.5\n    pa = landau_damping(rng, mesh, np, α, kx )\n    pm = ParticleMeshCoupling(pa, mesh)\n    energy = Float64[]\n    e = zeros(Float64, nx)\n    ρ = zeros(Float64, nx)\n    xmin = mesh.xmin\n    xmax = mesh.xmax\n    mat = compute_coeffs(pm, pa)\n    compute_rho!(ρ, mat, mesh, pa)\n    solve!(e, poisson, ρ)\n    for it in 1:nt+1       \n        update_positions!(pa, mesh, dt)\n        mat = compute_coeffs(pm, pa)\n        compute_rho!(ρ, mat, mesh, pa)\n        solve!(e, poisson, ρ)\n        update_velocities!(pa, e, mat, dt)\n        push!(energy, 0.5 * sum(e.^2) * mesh.dx) \n    end\n    energy\nend","category":"page"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"nt, dt = 1000, 0.01\nresults = main(nt, dt)\nt = collect(0:nt) .* dt\nplot( t, results, yaxis = :log )","category":"page"},{"location":"maxwell/#Maxwell-solver-using-Yee-scheme","page":"Maxwell solver","title":"Maxwell solver using Yee scheme","text":"","category":"section"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"L_xL_y domain dimensions and M,N are integers.","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"omega = sqrt(fracMpiL_x)^2+(fracNpiL_y)^2","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"B_z(xyt) =   - cos(M pi fracxL_x)  cos(N pi fracyL_y) cos(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"E_x(xyt) = fracc^2 N pi omega Ly cos(M pi fracxL_x) sin(N pi  fracyL_y) sin(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"E_y(xyt) = - fracc^2 M pi omega Lx sin (M pi fracxL_x) cos (N pi  fracyL_y) sin(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"using Plots","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"using ParticleInCell\n\ndimx, dimy = 1, 1\nnx, ny = 64, 64\nmd, nd = 2, 2  \ndt = 0.001\nnstep = 1 ÷ dt\nmesh = TwoDGrid( dimx, nx, dimy, ny )\nmaxwell = FDTD( mesh ) \nomega = sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)\n\nx = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])\ny = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose\n\nmaxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))\n    \nsurface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1))","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"Ex and Ey are set at t = 0.0\nBz is set at  t = -dt/2","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"function run(mesh, maxwell, nstep)\n\n    x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])\n    y = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose\n\n    maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))\n    \n    \n    @gif for istep = 1:nstep # Loop over time\n    \n        faraday!(maxwell, mesh, dt)     \n    \n        ampere_maxwell!(maxwell, mesh, dt) \n    \n        surface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1))\n\n    end every (nstep ÷ 100)\n    \n    \nend\n\nrun(mesh, maxwell, 2000)","category":"page"},{"location":"#ParticleInCell.jl","page":"Home","title":"ParticleInCell.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ParticleInCell.jl","category":"page"},{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Modules = [ParticleInCell]\nOrder   = [:type]","category":"page"},{"location":"types/#ParticleInCell.Maxwell1DFEM","page":"Types","title":"ParticleInCell.Maxwell1DFEM","text":"maxwell_solver = MaxwellFEM1D( mesh, degree )\n\n1D Maxwell spline finite element solver on a periodic grid\n\nLx                   : length of Periodic domain\ndelta_x              : cell size\nn_dofs               : number of cells (and grid points)\ns_deg_0              : spline degree 0-forms\ns_deg_1              : spline degree 1-forms\nmass_0               : coefficients of 0-form mass matrix\nmass_1               : coefficients of 1-form mass matrix\neig_mass0            : eigenvalues of circulant 0-form mass matrix\neig_mass1            : eigenvalues of circulant 1-form mass matrix\neig_weak_ampere      : eigenvalues of circulant update matrix for Ampere\neig_weak_poisson     : eigenvalues of circulant update matrix for Poisson\nplan_fw              : fft plan (forward)\nplan_bw              : fft plan (backward)\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.OneDGrid","page":"Types","title":"ParticleInCell.OneDGrid","text":"TwoDGrid( xmin, xmax, nx )\n\nSimple structure to store mesh data from 1 to 3 dimensions\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.PICPoisson2D","page":"Types","title":"ParticleInCell.PICPoisson2D","text":"PICPoisson2D( poisson, kernel )\n\nkernel : Kernel smoother taking care of charge deposition and field evaluation\npoisson : Poisson solver\nrho_dofs : Coefficients of expansion of rho (MPI global version)\nefield_dofs : Coefficients of expansion of electric field\nphi_dofs : Coefficients of expansion of potential\nrho2d : 2d version of rho_dofs to adjust to field solver format\nefield : 2d version of efield_dofs to adjust to field solver format\nphi2d : 2d version of phi_dofs to adjust to field solver format\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.ParticleGroup","page":"Types","title":"ParticleInCell.ParticleGroup","text":"ParticleGroup{D,V}( n_particles, charge, mass, q, weights)\n\nn_particles : number of particles \ncharge      : charge of the particle species\nmass        : mass of the particle species\nn_weights   : number of differents weights\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.ParticleMeshCoupling1D","page":"Types","title":"ParticleInCell.ParticleMeshCoupling1D","text":"ParticleMeshCoupling1D( mesh, no_particles, spline_degree, \n                      smoothing_type )\n\nKernel smoother with splines of arbitrary degree placed on a uniform mesh. Spline with index i starts at point i\n\ndelta_x : Value of grid spacing along both directions.\nn_grid : Array containing number ofpoints along each direction\nno_particles : Number of particles of underlying PIC method \nspline_degree : Degree of smoothing kernel spline\nn_span : Number of intervals where spline non zero (spline_degree + 1)\nscaling : Scaling factor depending on whether :galerkin or :collocation\nn_quad_points : Number of quadrature points\nspline_val: scratch data for spline evaluation\nspline_val_more : more scratch data for spline evaluation\nquad_x, quad_w : quadrature weights and points\n\nnote: Note\nOnly 1D version is implemented for now\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.ParticleMeshCoupling2D","page":"Types","title":"ParticleInCell.ParticleMeshCoupling2D","text":"ParticleMeshCoupling2D( pg, grid, degree)\n\nn_grid(2) : no. of spline coefficients\ndomain(2,2) : lower and upper bounds of the domain\nno_particles : no. of particles\ndegree : Degree of smoothing kernel spline\nsmoothing_type : Define if Galerkin or collocation smoothing for right scaling in accumulation routines \n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.Poisson2DPeriodic","page":"Types","title":"ParticleInCell.Poisson2DPeriodic","text":"Poisson2DPeriodic\n\nDerived type to solve the Poisson equation on 2d regular  cartesian mesh with periodic boundary conditions on both sides\n\nkx   : wave number in x\nky   : wave number in y\nk2   : k_x^2 + k_y^2\nnc_x : cells number in x\nnc_y : cells number in y\ndx   : x step size\ndy   : y step size\nrht  : fft(rho)\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.SplinePP","page":"Types","title":"ParticleInCell.SplinePP","text":"SplinePP( degree, ncells)\n\ndegree : degree of 1d spline\npoly_coeffs : poly_coeffs[i,j] coefficient of x^deg+1-j for ith B-spline function  size= (degree+1, degree+1)\npoly_coeffs_fp : poly_coeffs[i,j] coefficient of x^deg+1-j for ith B-spline function  size= (degree+1, degree+1)\nncells : number of gridcells\nscratch_b : scratch data for b_to_pp-converting\nscratch_p : scratch data for b_to_pp-converting\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.SplittingOperator","page":"Types","title":"ParticleInCell.SplittingOperator","text":"Operator splitting type for 2d2v Vlasov-Poisson\n\npic :: PIC poisson solver\npg :: Particle group\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.TwoDGrid","page":"Types","title":"ParticleInCell.TwoDGrid","text":"TwoDGrid( dimx, nx, dimy, ny)\n\nGenerate a cartesians mesh on rectangle dimxx dimy with nx x ny points\n\nnx : indices are in [1:nx]\nny : indices are in [1:ny]\ndimx = xmax - xmin\ndimy = ymax - ymin\nx, y : node positions\ndx, dy : step size\n\n\n\n\n\n","category":"type"}]
}
