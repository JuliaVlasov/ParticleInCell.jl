using FFTW

export Poisson2DPeriodic

"""
    Poisson2DPeriodic

Derived type to solve the Poisson equation on 2d regular 
cartesian mesh with periodic boundary conditions on both sides
- kx   : wave number in x
- ky   : wave number in y
- k2   : ``k_x^2 + k_y^2``
- nc_x : cells number in x
- nc_y : cells number in y
- dx   : x step size
- dy   : y step size
- rht  : fft(rho)
"""
struct Poisson2DPeriodic

   grid :: TwoDGrid
   kx   :: Array{Float64, 2}    
   ky   :: Array{Float64, 2}
   k2   :: Array{Float64, 2}
   rht  :: Array{ComplexF64, 2}

   function Poisson2DPeriodic( grid :: TwoDGrid )

       nc_x = grid.nx
       nc_y = grid.ny

       rht = zeros(ComplexF64,(div(nc_x,2)+1,nc_y))
    
       kx = zeros(nc_x÷2+1,nc_y)
       ky = zeros(nc_x÷2+1,nc_y)
       k2 = zeros(nc_x÷2+1,nc_y)

       kx0 = 2π/grid.dimx
       ky0 = 2π/grid.dimy
       
       for ik=1:nc_x÷2+1
          kx1 = (ik-1)*kx0
          for jk = 1:nc_y÷2
             kx[ik,jk] = kx1
             ky[ik,jk] = (jk-1)*ky0
          end
          for jk = nc_y÷2+1:nc_y     
             kx[ik,jk] = kx1
             ky[ik,jk] = (jk-1-nc_y)*ky0
          end
       end

       kx[1,1] = 1.0
       k2 .= kx .* kx .+ ky .* ky
       kx .= kx ./ k2
       ky .= ky ./ k2

       new( grid, kx, ky, k2, rht)

   end 

end

l2norm_squarred_2d_periodic(poisson, coefs_dofs) = sum(coefs_dofs.^2) * poisson.dx * poisson.dy

export solve!
    
"""
    solve!( poisson, phi, rho )

computes `phi` from `rho` 

```math
-\\Delta phi(x,y) = rho(x,y)
```

"""
function solve!( phi, poisson :: Poisson2DPeriodic, rho )
    
  poisson.rht .= rfft(rho)
  poisson.rht ./= poisson.k2
  phi .= irfft(poisson.rht, poisson.grid.nx)

end 

"""
    solve!( ex, ey, poisson, rho )

solves Poisson equation to compute electric fields

```math
E(x,y) = -\\nabla \\phi(x,y) \\\\
-\\Delta \\phi(x,y) = \\rho(x,y)
```

"""
function solve!( ex, ey, poisson :: Poisson2DPeriodic, rho )

  poisson.rht .= rfft(rho)

  poisson.rht[1,1] = 0.0
  poisson.rht .*=  -1im .* poisson.kx

  ex .= irfft(poisson.rht, poisson.grid.nx)  

  poisson.rht .= rfft(rho)

  poisson.rht[1,1] = 0.0
  poisson.rht .*=  -1im .* poisson.ky

  ey .= irfft(poisson.rht, poisson.grid.nx)  

end 


export PICPoisson2D

"""
    PICPoisson2D( poisson, kernel )

- kernel : Kernel smoother taking care of charge deposition and field evaluation
- poisson : Poisson solver
- rho_dofs : Coefficients of expansion of rho (MPI global version)
- efield_dofs : Coefficients of expansion of electric field
- phi_dofs : Coefficients of expansion of potential
- rho2d : 2d version of rho_dofs to adjust to field solver format
- efield : 2d version of efield_dofs to adjust to field solver format
- phi2d : 2d version of phi_dofs to adjust to field solver format
"""
struct PICPoisson2D
     
    ndofs :: Int
    kernel :: ParticleMeshCoupling2D 
    poisson :: Poisson2DPeriodic
    rho_dofs :: Vector{Float64}
    efield_dofs :: Vector{Vector{Float64}}
    phi_dofs :: Vector{Float64}
    rho2d :: Array{Float64, 2}
    efield :: Vector{Array{Float64, 2}}
    phi2d :: Array{Float64, 2}

    function PICPoisson2D( poisson :: Poisson2DPeriodic, 
                           kernel :: ParticleMeshCoupling2D)

        nx, ny = poisson.grid.nx, poisson.grid.ny
        ndofs = nx * ny
        
        rho_dofs = zeros(ndofs)
        rho_dofs_local = zeros(ndofs)
        rho_analyt_dofs = zeros(ndofs)
        efield_dofs = [zeros(ndofs), zeros(ndofs)]
        phi_dofs = zeros(ndofs)
        rho2d = zeros(nx, ny)
        efield = [zeros(nx, ny), zeros(nx, ny)]
        phi2d = zeros(nx, ny)

        new( ndofs, kernel, poisson, rho_dofs, 
             efield_dofs, phi_dofs, rho2d, efield, phi2d )

    end 

end

"""
    add_charge!(pic, position, marker_charge)

Add charge from one particle
- self : Pic Poisson solver object
- position : Position of the particle
- marker_charge : Particle weight times charge
"""
function add_charge!(pic :: PICPoisson2D, position, marker_charge)

    xp, yp = position

    add_charge!(pic.rho_dofs, pic.kernel, xp, yp, marker_charge)
       
end 
  
"""
    evaluate_rho!(pic, position)

Evaluate charge density at rho at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of rho at given position
"""
function evaluate_rho!(pic, position)

    evaluate!( kernel, position, pic.rho_dofs)

end

"""
    evaluate_phi!(pic, position)

Evaluate potential at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of phi at given position
"""
function evaluate_rho(pic, position)

    evaluate!( kernel, position, pic.phi_dofs)

end

"""
    evaluate_field_single_2d(pic, position, components, func_value)

Evaluate components of the electric field as one position
- self : PIC Poisson solver object
- components
- position : Position of the particle
- func_value
"""
evaluate_fields(pic, position) = evaluate_multiple( pic.kernel, position, pic.efield_dofs)


"""
    solve!( pic )

Solve for phi and fields
- poisson : Pic Poisson solver object
"""
function solve!(pic)

    solve_phi!(pic)
    solve_fields!(pic)
    
end

"""
    solve_phi!( pic )

Solve for potential
"""
function solve_phi!(pic)

    nx, ny = pic.poisson.grid.nx, pic.poisson.grid.ny
    pic.rho2d .= reshape(pic.rho_dofs, nx, ny)
    solve!(pic.phi2d, pic.poisson, pic.rho2d)
    pic.phi_dofs .= vec(pic.phi2d)

end 

"""
    solve_fields!( pic )

Solve efields from rho
"""
function solve_fields!(pic)

    nx, ny = pic.poisson.grid.nx, pic.poisson.grid.ny
    pic.rho2d .= reshape(pic.rho_dofs, nx, ny)
    solve!(pic.efield[1], pic.efield[2], pic.poisson, pic.rho2d)

    pic.efield_dofs[1] .= vec(pic.efield[1])
    pic.efield_dofs[2] .= vec(pic.efield[2])

end 


export SplittingOperator

"""
Operator splitting type for 2d2v Vlasov-Poisson
- pic :: PIC poisson solver
- pg :: Particle group
"""
struct SplittingOperator

    pic :: PICPoisson2D
    pg  :: ParticleGroup

end

"""
    operator_t(split, dt)
Push x 
- split :: time splitting object 
- dt   :: time step
"""
function operator_t!(split :: SplittingOperator, dt)

    # x_new = x_old + dt * v

    for i_part=1:split.pg.n_particles

        split.pg.array[1,i_part] += dt * split.pg.array[3,i_part] 
        split.pg.array[2,i_part] += dt * split.pg.array[4,i_part] 

    end

end
  
function operator_v!(split :: SplittingOperator, dt)

    # v_new = v_old + dt * q/m * E

    qm = split.pg.q_over_m

    for i_part=1:split.pg.n_particles

       # Evaluate efields at particle position

       xi = split.pg.array[1, i_part]
       yi = split.pg.array[2, i_part]

       ex, ey = evaluate_fields(split.pic, (xi, yi))

       vx_new = split.pg.array[3, i_part]
       vy_new = split.pg.array[4, i_part]

       vx_new += dt * qm * ex
       vy_new += dt * qm * ey

       split.pg.array[3,i_part] = vx_new
       split.pg.array[4,i_part] = vy_new

    end

end

export charge_deposition!

function charge_deposition!(split :: SplittingOperator)

    for i_part = 1:split.pg.n_particles
       xi = split.pg.array[1, i_part]
       yi = split.pg.array[2, i_part]
       wi = split.pg.array[5, i_part]
       add_charge!(split.pic, (xi, yi), wi)
    end

end

export solve_fields!

"""
Solve Poisson's equation for the electric field
"""
solve_fields!(split :: SplittingOperator) = solve_fields!( split.pic )


"""
Strang splitting
- split :: time splitting object 
- dt   :: time step
"""
function strang_splitting(split :: SplittingOperator, dt)

    operator_t!(split, 0.5dt)
    charge_deposition!(split)
    solve_fields!(split)
    operator_v!(split, dt)
    operator_t!(split, 0.5dt)

end
