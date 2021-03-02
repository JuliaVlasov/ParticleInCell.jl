using FFTW

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
    
       kx = zeros(nc_x/2+1,nc_y)
       ky = zeros(nc_x/2+1,nc_y)
       k2 = zeros(nc_x/2+1,nc_y)

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

export compute_phi_from_rho!
    
"""
    compute_phi_from_rho( poisson, phi, rho )

solves 

```math
-\\Delta phi(x,y) = rho(x,y)
```

"""
function compute_phi_from_rho!( phi, poisson, rho )
    
  poisson.rht .= rho
  fft!(poisson.rht)
  poisson.rht ./= poisson.k2
  ifft!(poisson.rht)
  phi .= real(poisson.rht)

end 

export compute_e_from_rho!

"""
    compute_e_from_rho!( ex, ey, poisson, rho )

solve Poisson equation to compute electric fields

```math
E(x,y) = -\\nabla \\phi(x,y) \\\\
-\\Delta \\phi(x,y) = \\rho(x,y)
```

"""
function compute_e_from_rho!( ex, ey, poisson, rho )

  poisson.rht = rho(1:nc_x,1:nc_y)
  fft!(poisson.rht)

  self.rht[1,1] = 0.0
  self.rht .*=  -1im .* poisson.kx

  ifft!(poisson.rht)  
  ex .= real(poisson.rht)

  poisson.rht = rho(1:nc_x,1:nc_y)
  fft!(poisson.rht)

  self.rht[1,1] = 0.0
  self.rht .*=  -1im .* poisson.ky

  ifft!(poisson.rht)  
  ey .= real(poisson.rht)

end 



"""
    PICPoisson2D

- kernel : Kernel smoother taking care of charge deposition and field evaluation
- poisson : Poisson solver
- rho_dofs : Coefficients of expansion of rho (MPI global version)
- efield_dofs : Coefficients of expansion of electric field
- phi_dofs : Coefficients of expansion of potential
- rho2d : 2d version of rho_dofs to adjust to field solver format
- efield1 : 2d version of efield_dofs(:,1) to adjust to field solver format
- efield2 : 2d version of efield_dofs(:,2) to adjust to field solver format
- phi2d : 2d version of phi_dofs to adjust to field solver format
"""
struct PICPoisson2D
     
    ndofs :: Int
    kernel :: ParticleMeshCoupling2D 
    poisson :: Poisson2DPeriodic
    rho_dofs :: Vector{Float64}
    efield_dofs :: Array{Float64, 2}
    phi_dofs :: Vector{Float64}
    rho2d :: Array{Float64, 2}
    efield1 :: Array{Float64, 2}
    efield2 :: Array{Float64, 2}
    phi2d :: Array{Float64, 2}

    function PICPoisson2D( poisson :: Poisson2DPeriodic, 
                           kernel :: ParticleMeshCoupling2D)

        nx, ny = poisson.grid.nx, poisson.grid.ny
        ndofs = nx * ny
        
        rho_dofs = zeros(ndofs)
        rho_dofs_local = zeros(ndofs)
        rho_analyt_dofs = zeros(ndofs)
        efield_dofs = zeros(ndofs, 2)
        phi_dofs = zeros(ndofs)
        rho2d = zeros(nx, ny)
        efield1 = zeros(nx, ny)
        efield2 = zeros(nx, ny)
        phi2d = zeros(nx, ny)

        new( ndofs, kernel, poisson, rho_dofs, 
             efield_dofs, phi_dofs, rho2d, efield1, efield2, phi2d )

    end 

end

"""
    add_charge_single_2d!(pic, position, marker_charge)

Add charge from one particle
- self : Pic Poisson solver object
- position : Position of the particle
- marker_charge : Particle weight times charge
"""
function add_charge_single_2d(pic, position, marker_charge)

    add_charge!(pic.rho_dofs_local, kernel, position, marker_charge)
       
end 
  
"""
    evaluate_rho_single_2d(self, position, func_value)

Evaluate charge density at rho at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of rho at given position
"""
function evaluate_rho_single_2d(pic, position)

    evaluate( kernel, position, pic.rho_dofs)

end

"""
    evaluate_phi_single_2d(self, position, func_value)

Evaluate potential at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of phi at given position
"""
function evaluate_phi_single_2d(pic, position)

    evaluate( kernel, position, pic.phi_dofs)

end

"""
    evaluate_field_single_2d(pic, position, components, func_value)

Evaluate components of the electric field as one position
- self : PIC Poisson solver object
- components
- position : Position of the particle
- func_value
"""
function evaluate_field_single_2d(pic, position, components)

    evaluate_multiple( kernel, position, components, pic.efield_dofs)

end


"""
    solve!( pic )

Solve for phi and fields
- poisson : Pic Poisson solver object
"""
function solve!(pic)

    solve_phi(pic)
    solve_fields(pic)
    
end

"""
    solve_phi!( pic )

Solve for potential
"""
function solve_phi!(pic)

    nx, ny = pic.poisson.grid.nx, pic.poisson.grid.ny
    self.rho2d .= reshape(pic.rho_dofs, nx, ny)
    compute_phi_from_rho!(pic.phi2d, pic.poisson, pic.rho2d)
    pic.phi_dofs .= vec(pic.phi2d)

end 

"""
    solve_fields!( pic )

Solve efields from rho
"""
function solve_fields!(pic)

    nx, ny = poisson.grid.nx, poisson.grid.ny
    pic.rho2d = reshape(pic.rho_dofs, nx, ny)
    compute_e_from_rho!(pic.efield1, pic.efield2, pic.poisson, pic.rho2d)

    pic.efield_dofs[1] .= vec(pic.efield1)
    pic.efield_dofs[2] .= vec(pic.efield2)

end 
