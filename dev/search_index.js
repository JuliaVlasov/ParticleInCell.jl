var documenterSearchIndex = {"docs":
[{"location":"vlasov-maxwell/#Vlasov-Maxwell","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell","text":"","category":"section"},{"location":"vlasov-maxwell/#Landau-damping","page":"Vlasov-Maxwell 2D","title":"Landau damping","text":"","category":"section"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"using Plots, LinearAlgebra","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"using ParticleInCell\n\nconst nx = 128 \nconst ny = 16   \n\nalpha = 0.1\nkx = 0.5\nky = 0.\ndimx = 2*pi/kx\ndimy = 1  \npoids = dimx * dimy \n\nmesh = TwoDGrid( dimx, nx, dimy, ny)\nfdtd = FDTD(mesh)\n\n\ntime  = 0\n\nfor i=1:nx, j=1:ny+1\n    fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)\nend\nsurface(fdtd.ex )","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"nbpart = 100*nx*ny\nparticles = ParticleGroup{2,2}( nbpart, charge=1.0, mass=1.0, n_weights=1)\nsampler = LandauDamping( alpha, kx )\nsample!( particles, mesh, sampler)\n\np = plot(layout=4)\nhistogram!(p[1], particles.array[1,:], normalize=true, label=\"x\")\nhistogram!(p[2], particles.array[2,:], normalize=true, label=\"y\")\nhistogram!(p[3], particles.array[3,:], normalize=true, label=\"vx\")\nhistogram!(p[4], particles.array[4,:], normalize=true, label=\"vy\")","category":"page"},{"location":"vlasov-maxwell/","page":"Vlasov-Maxwell 2D","title":"Vlasov-Maxwell 2D","text":"function run( nstep; npm = 100 )\n    \n    dt = 0.01\n    alpha = 0.1\n    kx   = 0.5\n    dimx = 2*pi/kx\n    dimy = 1  \n    nx   = 128  # nombre de pts suivant x\n    ny   = 16   # nombre de pts suivant y\n    mesh = TwoDGrid( dimx, nx, dimy, ny)\n    dx, dy = mesh.dx, mesh.dy\n\n    ex = zeros(nx+1, ny+1)\n    ey = zeros(nx+1, ny+1)\n    bz = zeros(nx+1, ny+1)\n    jx = zeros(nx+1, ny+1)\n    jy = zeros(nx+1, ny+1)\n\n    nbpart = npm*nx*ny\n    println( \" nbpart = $nbpart \")\n\n    particles = ParticleGroup{2,2}( nbpart, charge=1.0, mass=1.0, n_weights=1)\n    sampler = LandauDamping( alpha, kx )\n    sample!( particles, mesh, sampler)\n\n    fdtd = FDTD(mesh)\n    for i=1:nx, j=1:ny+1\n        fdtd.ex[i,j] = alpha/kx * sin(kx*(mesh.x[i]+mesh.x[i+1])/2)\n    end\n    time = 0\n    energy = Float64[compute_energy(fdtd, mesh)]\n    t = Float64[time]\n\n    kernel = CloudInCell()\n    \n    for istep in 1:nstep\n    \n       if istep > 1\n           faraday!( fdtd, mesh, 0.5dt ) \n       end\n       update_fields!(ex, ey, bz, mesh, fdtd)\n       push_v!( particles, kernel, mesh, ex, ey, bz, dt )\n       push_x!( particles, mesh, 0.5dt) \n       compute_current!( jx, jy, mesh, kernel, particles)\n       push_x!( particles, mesh, 0.5dt) \n       faraday!(fdtd, mesh, 0.5dt)\n       ampere_maxwell!(fdtd, mesh, jx, jy, dt)\n       time = time + dt\n       push!(t, time)\n       push!(energy, compute_energy(fdtd, mesh))\n    \n    end\n   \n    t, energy\n    \nend\n\nnstep = 1000\nt, energy = run(nstep)\nplot(t, energy)","category":"page"},{"location":"functions/#Functions","page":"Functions","title":"Functions","text":"","category":"section"},{"location":"functions/","page":"Functions","title":"Functions","text":"Modules = [ParticleInCell]\nOrder   = [:function]","category":"page"},{"location":"functions/#ParticleInCell.add_charge!-Tuple{TwoDPoissonPIC,Any,Any}","page":"Functions","title":"ParticleInCell.add_charge!","text":"add_charge!(pic, position, marker_charge)\n\nAdd charge from one particle\n\nself : Pic Poisson solver object\nposition : Position of the particle\nmarker_charge : Particle weight times charge\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.compute_field_energy-Tuple{TwoDPoissonPIC,Int64}","page":"Functions","title":"ParticleInCell.compute_field_energy","text":"compute_field_energy\n\nCompute the squared l2 norm of electric field\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_rho!-Tuple{Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_rho!","text":"evaluate_rho!(pic, position)\n\nEvaluate charge density at rho at one position\n\nself : Pic Poisson solver object\nposition : Position of the particle\nfunc_value : Value of rho at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.evaluate_rho-Tuple{Any,Any}","page":"Functions","title":"ParticleInCell.evaluate_rho","text":"evaluate_phi!(pic, position)\n\nEvaluate potential at one position\n\nself : Pic Poisson solver object\nposition : Position of the particle\nfunc_value : Value of phi at given position\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.operator_t!-Tuple{SplittingOperator,Any}","page":"Functions","title":"ParticleInCell.operator_t!","text":"operator_t(split, dt)\n\nPush x \n\nsplit :: time splitting object \ndt   :: time step\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.sample!-Tuple{ParticleGroup{1,1},OneDGrid,LandauDamping}","page":"Functions","title":"ParticleInCell.sample!","text":"sample!(d, pg)\n\nSampling from a probability distribution to initialize a Landau damping in 1D1V space.\n\nf_0(xvt) = fracn_02π v_th^2 ( 1 + alpha cos(k_x x)) exp( - fracv^22 v_th^2)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.sample!-Tuple{ParticleGroup{1,2},LandauDamping}","page":"Functions","title":"ParticleInCell.sample!","text":"sample!(d, pg)\n\nSampling from a probability distribution to initialize a Landau damping in 1D2V space.\n\nf_0(xvt) = fracn_02π v_th^2 ( 1 + alpha cos(k_x x))\n exp( - fracv_x^2+v_y^22 v_th^2)\n\nThe newton function solves the equation P(x)-r=0 with Newton’s method\n\nx^n+1 = x^n  (P(x)-(2pi r  k)f(x) \n\nwith \n\nP(x) = int_0^x (1 + alpha cos(k_x y)) dy = x + fracalphak_x sin(k_x x)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.sample!-Tuple{ParticleGroup{2,2},TwoDGrid,LandauDamping}","page":"Functions","title":"ParticleInCell.sample!","text":"sample!(d, pg)\n\nSampling from a probability distribution to initialize a Landau damping in 2D2V space.\n\nf_0(xvt) = fracn_02π v_th^2 ( 1 + alpha cos(k_x x))\n exp( - fracv_x^2+v_y^22 v_th^2)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any,Any,TwoDPoissonPeriodic,Any}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( ex, ey, poisson, rho )\n\nsolves Poisson equation to compute electric fields\n\nE(xy) = -nabla phi(xy) \n-Delta phi(xy) = rho(xy)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{Any,TwoDPoissonPeriodic,Any}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( poisson, phi, rho )\n\ncomputes phi from rho \n\n-Delta phi(xy) = rho(xy)\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve!-Tuple{TwoDPoissonPIC}","page":"Functions","title":"ParticleInCell.solve!","text":"solve!( pic )\n\nSolve for phi and fields\n\npoisson : Pic Poisson solver object\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_fields!-Tuple{SplittingOperator}","page":"Functions","title":"ParticleInCell.solve_fields!","text":"Solve Poisson's equation for the electric field\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_fields!-Tuple{TwoDPoissonPIC}","page":"Functions","title":"ParticleInCell.solve_fields!","text":"solve_fields!( pic )\n\nSolve efields from rho\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.solve_phi!-Tuple{TwoDPoissonPIC}","page":"Functions","title":"ParticleInCell.solve_phi!","text":"solve_phi!( pic )\n\nSolve for potential\n\n\n\n\n\n","category":"method"},{"location":"functions/#ParticleInCell.strang_splitting!-Tuple{SplittingOperator,Any}","page":"Functions","title":"ParticleInCell.strang_splitting!","text":"Strang splitting\n\nsplit :: time splitting object \ndt   :: time step\n\n\n\n\n\n","category":"method"},{"location":"tsi/#Two-stream-instability","page":"Two-stream instability","title":"Two-stream instability","text":"","category":"section"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"using Plots\nusing Random","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"using ParticleInCell\n\nconst dt = 0.005     # Time step\nconst nt = 10000     # Number of time steps\nconst L  = 20π       #  Domain size \nconst nx = 320       # Number of grid cells\nconst np = nx * 20   # Number of particles\n\n\nmesh = OneDGrid( 0, 20π, nx)\nrng = MersenneTwister(42)\npoisson = OneDPoisson( mesh )\nparticles = tsi(rng, mesh, np )\npm = ParticleMeshCoupling(particles, mesh)","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"function main()\n\n    mesh = OneDGrid( 0, 20π, nx)\n    poisson = OneDPoisson( mesh )\n    rng = MersenneTwister(42)\n    pa = tsi(rng, mesh, np )\n    pm = ParticleMeshCoupling(pa, mesh)\n    energy = Float64[]\n    e = zeros(Float64, nx)\n    ρ = zeros(Float64, nx)\n    xmin = mesh.xmin\n    xmax = mesh.xmax\n    \n    for it in 1:nt+1\n        \n        update_positions!(pa, mesh, dt)\n        mat = compute_coeffs(pm, pa)\n        compute_rho!(ρ, mat, mesh, pa)\n        solve!(e, poisson, ρ)\n        update_velocities!(pa, e, mat, dt)\n        push!(energy, 0.5 * sum(e.^2) * mesh.dx) \n\n    end\n\n    energy\n\nend","category":"page"},{"location":"tsi/","page":"Two-stream instability","title":"Two-stream instability","text":"results = main()\nt = (0:nt) .* dt\nplot( t, results, yaxis=:log)","category":"page"},{"location":"contents/#Contents","page":"Contents","title":"Contents","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"contents/#Index","page":"Contents","title":"Index","text":"","category":"section"},{"location":"contents/","page":"Contents","title":"Contents","text":"","category":"page"},{"location":"landau_damping/#Landau-damping","page":"Landau damping","title":"Landau damping","text":"","category":"section"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"using ParticleInCell\nusing Plots\nusing Random\nusing GEMPIC","category":"page"},{"location":"landau_damping/#Test-particle-initialization","page":"Landau damping","title":"Test particle initialization","text":"","category":"section"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"dt = 0.1\nnsteps = 100\nalpha = 0.1\nkx = 0.5\n\nnx = 128\nxmin, xmax = 0.0, 2π / kx\n\nn_particles = 100000\ndegree_smoother = 3\n\nmesh = OneDGrid( xmin, xmax, nx)\n\nparticles = ParticleGroup{1,1}( n_particles, charge=1.0, mass=1.0, n_weights=1)\n\nsampler = LandauDamping( alpha, kx )\n\nParticleInCell.sample!( particles, mesh, sampler)\n\nparticles.array[3,:] .= (xmax - xmin) ./ n_particles;","category":"page"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"p = plot(layout=2)\nhistogram!( p[1], particles.array[1,:], normalized=true)\nhistogram!( p[2], particles.array[2,:], normalized=true)","category":"page"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"poisson = OneDPoisson( mesh )\nkernel = ParticleMeshCoupling1D( particles, mesh, degree_smoother, :collocation)\n\nex = zeros(nx)\nrho = zeros(nx)\n\nfor i_part = 1:particles.n_particles\n    xi = particles.array[1, i_part]\n    wi = particles.array[3, i_part]\n    GEMPIC.add_charge!(rho, kernel, xi, wi)\nend\n\nsolve!(ex, poisson, rho)\nx = LinRange( xmin, xmax, nx+1)[1:end-1]\np = plot(layout=(2))\nplot!(p[1], x, ex)\nplot!(p[1], x, alpha/kx * sin.(kx * x))\nplot!(p[2], x, rho)\nplot!(p[2], x, alpha * cos.(kx * x))","category":"page"},{"location":"landau_damping/","page":"Landau damping","title":"Landau damping","text":"\nproblem = OneDPoissonPIC( poisson, kernel )\n\ndt = 0.1\nnsteps = 100\nalpha = 0.1\nkx = 0.5\n\npropagator = SplittingOperator( problem, particles ) \n\nenergy = Float64[]\n\nfor j=1:nsteps\n\n    ParticleInCell.operator_t!(propagator, 0.5dt)\n    ParticleInCell.charge_deposition!(propagator)\n    ParticleInCell.solve_fields!(propagator)\n    ParticleInCell.operator_v!(propagator, dt)\n    ParticleInCell.operator_t!(propagator, 0.5dt)\n\n    push!(energy, compute_field_energy(problem, 1))\n          \nend\n\nt = collect(0:nsteps) .* dt\nplot(log.(energy))\n","category":"page"},{"location":"maxwell/#Maxwell-solver-using-Yee-scheme","page":"Maxwell solver","title":"Maxwell solver using Yee scheme","text":"","category":"section"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"L_xL_y domain dimensions and M,N are integers.","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"omega = sqrt(fracMpiL_x)^2+(fracNpiL_y)^2","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"B_z(xyt) =   - cos(M pi fracxL_x)  cos(N pi fracyL_y) cos(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"E_x(xyt) = fracc^2 N pi omega Ly cos(M pi fracxL_x) sin(N pi  fracyL_y) sin(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"E_y(xyt) = - fracc^2 M pi omega Lx sin (M pi fracxL_x) cos (N pi  fracyL_y) sin(omega t)","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"using Plots","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"using ParticleInCell\n\ndimx, dimy = 1, 1\nnx, ny = 64, 64\nmd, nd = 2, 2  \ndt = 0.001\nnstep = 1 ÷ dt\nmesh = TwoDGrid( dimx, nx, dimy, ny )\nmaxwell = FDTD( mesh ) \nomega = sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)\n\nx = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])\ny = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose\n\nmaxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))\n    \nsurface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1))","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"Ex and Ey are set at t = 0.0\nBz is set at  t = -dt/2","category":"page"},{"location":"maxwell/","page":"Maxwell solver","title":"Maxwell solver","text":"function run(mesh, maxwell, nstep)\n\n    x = 0.5 .* (mesh.x[1:end-1] .+ mesh.x[2:end])\n    y = 0.5 .* (mesh.y[1:end-1] .+ mesh.y[2:end]) |> transpose\n\n    maxwell.bz .= - cos.(md*pi*x) .* cos.(nd*pi*y) .* cos(omega*(-0.5*dt))\n    \n    nx, ny = mesh.nx, mesh.ny\n    jx = zeros(nx+1, ny+1)\n    jy = zeros(nx+1, ny+1)\n    \n    @gif for istep = 1:nstep # Loop over time\n    \n        faraday!(maxwell, mesh, dt)     \n    \n        ampere_maxwell!(maxwell, mesh, jx, jy, dt) \n    \n        surface(maxwell.bz, aspect_ratio=:equal, zlims=(-1,1), clim=(-1,1))\n\n    end every (nstep ÷ 100)\n    \n    \nend\n\nrun(mesh, maxwell, 2000)","category":"page"},{"location":"#ParticleInCell.jl","page":"Home","title":"ParticleInCell.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for ParticleInCell.jl","category":"page"},{"location":"vlasov-poisson/#Vlasov-Poisson","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson","text":"","category":"section"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"Simulation of 2d2v Vlasov-Poisson with simple PIC method, periodic boundary conditions and Landau initial values along x1 only.","category":"page"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"using Plots\nusing GEMPIC\nusing ParticleInCell\n\ndt = 0.1\nnsteps = 100\nalpha = 0.1\nkx = 0.5\n\nnx = 128\nny = 16\nxmin, xmax = 0.0, 4π\nymin, ymax = 0.0, 1.0\n\nn_particles = 100000\ndegree_smoother = 3\n\nmesh = TwoDGrid( xmin, xmax, nx, ymin, ymax, ny)\n\nparticles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)\n\nsampler = LandauDamping( alpha, kx )\n\nParticleInCell.sample!( particles, mesh, sampler)\n\nparticles.array[5,:]  .= (xmax - xmin) * (ymax - ymin) ./ n_particles;","category":"page"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"p = plot(layout=2)\nhistogram!( p[1], particles.array[1,:], normalized=true)\nhistogram!( p[2], particles.array[2,:], normalized=true)","category":"page"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"p = plot(layout=2)\n\nhistogram!( p[1], particles.array[3,:], normalized=true)\nxlims!(-6,6)\nhistogram!( p[2], particles.array[4,:], normalized=true)\nxlims!(-6,6)","category":"page"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"poisson = TwoDPoissonPeriodic( mesh )\nkernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)\n\nex = zeros(nx, ny)\ney = zeros(nx, ny)\nrho_dofs = zeros(nx*ny)\n\nfor i_part = 1:particles.n_particles\n    xi = particles.array[1, i_part]\n    yi = particles.array[2, i_part]\n    wi = particles.array[5, i_part]\n    GEMPIC.add_charge!(rho_dofs, kernel, xi, yi, wi)\nend\nrho = reshape(rho_dofs, nx, ny )\nsolve!(ex, ey, poisson, rho)\np = plot(layout=(2))\nsurface!(p[1], ex)\nsurface!(p[2],rho)","category":"page"},{"location":"vlasov-poisson/","page":"Vlasov-Poisson 2D","title":"Vlasov-Poisson 2D","text":"\npoisson = TwoDPoissonPeriodic( mesh )\n\nkernel = ParticleMeshCoupling2D( particles, mesh, degree_smoother, :collocation)\n\nproblem = TwoDPoissonPIC(  poisson, kernel )\n\ndt = 0.1\nnsteps = 100\nalpha = 0.1\nkx = 0.5\n\nparticles = ParticleGroup{2,2}( n_particles, charge=1.0, mass=1.0, n_weights=1)\n\nsampler = LandauDamping( alpha, kx )\n\nParticleInCell.sample!( particles, mesh, sampler)\n\nparticles.array[2,:] .*= ( ymax - ymin)\nparticles.array[5,:]  .= (xmax - xmin) * (ymax - ymin) / n_particles;\n\npropagator = SplittingOperator( problem, particles ) \n\nenergy = Float64[]\n\nfor j=1:100\n\n    ParticleInCell.operator_t!(propagator, 0.5dt)\n    ParticleInCell.charge_deposition!(propagator)\n    ParticleInCell.solve_fields!(propagator)\n    ParticleInCell.operator_v!(propagator, dt)\n    ParticleInCell.operator_t!(propagator, 0.5dt)\n\n    push!(energy, compute_field_energy(problem, 1))\n          \nend\n\nplot(log.(energy))\n\n","category":"page"},{"location":"types/#Types","page":"Types","title":"Types","text":"","category":"section"},{"location":"types/","page":"Types","title":"Types","text":"Modules = [ParticleInCell]\nOrder   = [:type]","category":"page"},{"location":"types/#ParticleInCell.LandauDamping","page":"Types","title":"ParticleInCell.LandauDamping","text":"Landau( α, kx)\n\nTest structure to initialize a particles distribtion for Landau damping test case in 1D1V and 1D2V\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.SplittingOperator","page":"Types","title":"ParticleInCell.SplittingOperator","text":"Operator splitting type for 2d2v Vlasov-Poisson\n\npic :: PIC poisson solver\npg :: Particle group\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.TwoDPoissonPIC","page":"Types","title":"ParticleInCell.TwoDPoissonPIC","text":"TwoDPoissonPIC( poisson, kernel )\n\nkernel : Kernel smoother taking care of charge deposition and field evaluation\npoisson : Poisson solver\nrho_dofs : Coefficients of expansion of rho (MPI global version)\nefield_dofs : Coefficients of expansion of electric field\nphi_dofs : Coefficients of expansion of potential\nrho2d : 2d version of rho_dofs to adjust to field solver format\nefield : 2d version of efield_dofs to adjust to field solver format\nphi2d : 2d version of phi_dofs to adjust to field solver format\n\n\n\n\n\n","category":"type"},{"location":"types/#ParticleInCell.TwoDPoissonPeriodic","page":"Types","title":"ParticleInCell.TwoDPoissonPeriodic","text":"TwoDPoissonPeriodic\n\nDerived type to solve the Poisson equation on 2d regular  cartesian mesh with periodic boundary conditions on both sides\n\nkx   : wave number in x\nky   : wave number in y\nk2   : k_x^2 + k_y^2\nnc_x : cells number in x\nnc_y : cells number in y\ndx   : x step size\ndy   : y step size\nrht  : fft(rho)\n\n\n\n\n\n","category":"type"}]
}
