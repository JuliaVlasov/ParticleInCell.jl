using FFTW

```
    Poisson2DPeriodic

derived type to solve the Poisson equation on 2d regular 
cartesian mesh with periodic boundary conditions on both sides
- kx       : wave number in x
- ky       : wave number in y
- k2       : ``k_x^2 + k_y^2 ``
- nc_x     : cells number in x
- nc_y     : cells number in y
- dx       : x step size
- dy       : y step size
- rht(:,:) : fft(rho)
- exy(:,:) : fft(ex and ey)
- fw       : forward fft plan
- bw       : backward fft plan
- p_rho    : C array pointer
- p_exy    : C array pointer
- p_tmp    : C array pointer
- tmp(:,:)
```
struct Poisson2DPeriodic

   kx   :: Array{Float64, 2}    
   ky   :: Array{Float64, 2}
   k2   :: Array{Float64, 2}
   rht  :: Array{ComplexF64, 2}
   exy  :: Array{ComplexF64, 2}
   fw   :: FFTWPlan    
   bw   :: FFTWPlan    
   tmp  :: Array{Float64, 2}

   function Poisson2DPeriodic( grid :: TwoDGrid )

       rht(1:nc_x/2+1,1:nc_y)
       exy(1:nc_x/2+1,1:nc_y)
       tmp(1:nc_x,1:nc_y)
       fw = fft_plan( tmp, rht)
       bw = fft_plan( rht, tmp )
        
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

   end 

end

l2norm_squarred_2d_periodic(poisson, coefs_dofs) = sum(coefs_dofs.^2) * poisson.dx * poisson.dy
    
```
solves $$ -\Delta phi(x,y) = rho (x,y) $$
```
function compute_phi_from_rho( poisson, phi, rho )
    
    solve( poisson, phi, rho)
    
end 

```
solve Poisson equation to compute electric fields

```math
E(x,y) = -\\nabla \\phi(x,y) \\\\
-\D\elta \\phi(x,y) = \\rho(x,y)
```

```
function compute_e_from_rho( poisson, e1, e2, rho )

    solve( poisson%solver, E1, E2, rho)
      
end 


function solve(self, phi, rho)

  nc_x = self.grid.nx
  nc_y = self.grid.ny

  self.tmp .= rho
  tmp .= fw * self.rht

  rht .= rht ./ self.k2

  self.rht .= bw * self.tmp

end

function solve(self,e_x,e_y,rho,nrj)

  self.tmp = rho(1:nc_x,1:nc_y)
  call sll_s_fft_exec_r2c_2d(self.fw, self%tmp, self%rht)

  self.exy(1,1) = (0.0_f64,0.0_f64)
  self.exy = -cmplx(0.0_f64,self%kx,kind=f64)*self%rht
  call sll_s_fft_exec_c2r_2d(self.bw, self%exy, self%tmp)
  e_x(1:nc_x,1:nc_y) = self.tmp / real(nc_x*nc_y, f64)

  self.exy(1,1) = (0.0_f64,0.0_f64)
  self.exy = -cmplx(0.0_f64,self%ky,kind=f64)*self%rht
  call sll_s_fft_exec_c2r_2d(self.bw, self%exy, self%tmp)

  e_y(1:nc_x,1:nc_y) = self.tmp / real(nc_x*nc_y, f64)

  !Node centered case
  if (size(e_x,1) == nc_x+1) e_x(nc_x+1,:) = e_x(1,:)
  if (size(e_x,2) == nc_y+1) e_x(:,nc_y+1) = e_x(:,1)
  if (size(e_y,1) == nc_x+1) e_y(nc_x+1,:) = e_y(1,:)
  if (size(e_y,2) == nc_y+1) e_y(:,nc_y+1) = e_y(:,1)

  if (present(nrj)) then 
     dx = self.dx
     dy = self.dy
     nrj=sum(e_x(1:nc_x,1:nc_y)*e_x(1:nc_x,1:nc_y) &
       +e_y(1:nc_x,1:nc_y)*e_y(1:nc_x,1:nc_y))*dx*dy
  end

end function solve_e_fields_poisson_2d_periodic_fft

!> Delete the Poisson object
function delete_poisson_2d_periodic_fft(self)

  type(sll_t_poisson_2d_periodic_fft) :: self

  call sll_s_fft_free(self.fw)
  call sll_s_fft_free(self.bw)

end function delete_poisson_2d_periodic_fft

end module sll_m_poisson_2d_periodic
```
    PICPoisson2D

- kernel : Kernel smoother taking care of charge deposition and field evaluation
- poisson : Poisson solver
- rho_dofs : Coefficients of expansion of rho (MPI global version)
- rho_dofs_local : Coefficients of expansion of rho (MPI local version)
- rho_analyt_dofs : Analytic contribution to the coefficients of expansion of rho (for delta f)
- efield_dofs : Coefficients of expansion of electric field
- phi_dofs : Coefficients of expansion of potential
- rho2d : 2d version of rho_dofs to adjust to field solver format
- efield1 : 2d version of efield_dofs(:,1) to adjust to field solver format
- efield2(:,:) : 2d version of efield_dofs(:,2) to adjust to field solver format
- phi2d(:,:) : 2d version of phi_dofs to adjust to field solver format
- rho_collected : Flag to indicate if charge deposition has been finished
```
struct PICPoisson2D
     
    no_gridpts :: Int
    no_dofs :: Int
    kernel :: ParticleMeshCOupling2D 
    poisson :: Poisson
    rho_dofs :: Vector{Float64}
    rho_dofs_local :: Vector{Float64}
    rho_analyt_dofs :: Vector{Float64}
    efield_dofs :: Array{Float64, 2}
    phi_dofs :: Vector{Float64}
    rho2d :; Array{Float64, 2}
    efield1 :; Array{Float64, 2}
    efield2 :: Array{Float64, 2}
    phi2d :: Array{Float64, 2}
    rho_collected :: Bool


    function PICPoisson2D( no_gridpts, poisson, kernel)

        dim = 1
        no_gridpts = no_gridpts
        no_dofs = prod(no_gridpts)
        rho_collected = false
        
        rho_dofs = zeros(no_dofs)
        rho_dofs_local = zeros(no_dofs)
        rho_analyt_dofs = zeros(no_dofs)
        efield_dofs = zeros(no_dofs, 2)
        phi_dofs = zeros(no_dofs)
        rho2d = zeros(no_gridpts(1), no_gridpts(2))
        efield1 = zeros(no_gridpts(1), no_gridpts(2))
        efield2 = zeros(no_gridpts(1), no_gridpts(2))
        phi2d = zeros(no_gridpts(1), no_gridpts(2))

        new( no_gridpts, no_dofs, kernel, poisson, rho_dofs, rho_dofs_local,
             rho_analyt_dofs, efield_dofs, phi_dofs, rho2d, efield1, efield2,
             phi2d, rho_collected)

    end 

end

"""
    add_charge_single_2d(pic, position, marker_charge)
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
function evaluate_rho_single_2d(self, position, func_value)

    if !self.rho_collected
       rho_collected = true
       fill!(self.rho_dofs, 0.0)
       self.rho_dofs .= reduce(+, self.rho_dofs_local)
    end

    func_value = evaluate( kernel position, sel.%rho_dofs)

end

```
    evaluate_phi_single_2d(self, position, func_value)

Evaluate potential \a phi at one position
- self : Pic Poisson solver object
- position : Position of the particle
- func_value : Value of phi at given position
```
function evaluate_phi_single_2d(self, position, func_value)

    func_value = evaluate( kernel, position, self.phi_dofs)

end

```
    evaluate_field_single_2d(self, position, components, func_value)
Evaluate components \a components of the electric field as one position
- self : Pic Poisson solver object
- components
- position : Position of the particle
- func_value
```
function evaluate_field_single_2d(self, position, components, func_value)

    func_value = evaluate_multiple( kernel, position, components, self.efield_dofs)

end


```
    solve( poisson )

Solve for phi and fields
- poisson : Pic Poisson solver object
```
function solve(self)

    solve_phi(self)
    solve_fields(self)
    
end

```
Solve for potential
```
function solve_phi_2d(self)

    if (self.rho_collected .EQV. .FALSE.) then
       self.rho_collected = .TRUE.
       self.rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self.rho_dofs_local, &
            self.no_dofs, MPI_SUM, self%rho_dofs)
    end
    self.rho2d = reshape(self%rho_dofs, self%no_gridpts)
    call self.poisson%compute_phi_from_rho(self%phi2d, self%rho2d)
    self.phi_dofs = reshape(self%phi2d, [self%no_dofs])

end 

```
Solve efields from rho
```
function solve_fields_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    if (self.rho_collected .EQV. .FALSE.) then
       self.rho_collected = .TRUE.
       self.rho_dofs = 0.0_f64
       call sll_o_collective_allreduce( sll_v_world_collective, self.rho_dofs_local, &
            self.no_dofs, MPI_SUM, self%rho_dofs)
    end
    self.rho2d = reshape(self%rho_dofs, self%no_gridpts)
    call self.poisson%compute_E_from_rho(self%efield1, self%efield2, self%rho2d)
    self.efield_dofs(:,1) = reshape(self%efield1, [self%no_dofs])
    self.efield_dofs(:,2) = reshape(self%efield2, [self%no_dofs])

end 

```
Reset charge to zero
```
function reset_2d(self)
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object

    self.rho_dofs_local = 0.0_f64
    self.rho_dofs = 0.0_f64
    self.rho_collected = .FALSE.

end 

```
Add analytic charge (set by \a set_analytic_charge ) to the accumulated charge
```
  function add_analytic_charge_2d(self, factor_present, factor_analytic)   
    class(sll_t_pic_poisson_2d), intent( inout ) :: self !< Pic Poisson solver object
    sll_real64, intent( in ) :: factor_present !< Factor to multiply accumulated charge with
    sll_real64, intent( in ) :: factor_analytic !< Factor to multiply added analytic charge with

    self.rho_dofs = factor_present * self%rho_dofs + &
         factor_analytic * self.rho_analyt_dofs

  end 

  !> Set analytic charge defined by a function \a func obeying the interface \a sll_i_function_of_position
  function set_analytic_charge_2d(self, func)
    class( sll_t_pic_poisson_2d ), intent( inout )    :: self !< PIC Poisson solver object.
    procedure(sll_i_function_of_position)                :: func !< Function to be projected.

    call self.poisson%compute_rhs_from_function(func, self%rho_analyt_dofs)

  end 

  !> Compute the squared l2 norm of component \a component of the field
  function compute_field_energy_2d(self, component) result(energy)
    class (sll_t_pic_poisson_2d), intent( in ) :: self
    sll_int32,                    intent( in ) :: component !< Component of the electric field for which the energy should be computed
    sll_real64                                 :: energy !< L2 norm squarred of 


    if (component == 1) then
       energy = self.poisson%l2norm_squared(self%efield1)
    elseif (component == 2) then
       energy = self.poisson%l2norm_squared(self%efield2)
    end

  end 

