export ParticleMeshCoupling2D

"""
    ParticleMeshCoupling2D( pg, grid, spline_degree)
- n_grid(2) : no. of spline coefficients
- domain(2,2) : lower and upper bounds of the domain
- no_particles : no. of particles
- spline_degree : Degree of smoothing kernel spline
- smoothing_type : Define if Galerkin or collocation smoothing for right scaling in accumulation routines 
"""
struct ParticleMeshCoupling2D

    grid :: TwoDGrid
    npart :: Int
    spline :: Vector{SplinePP}

    function ParticleMeshCoupling2D( pg :: ParticleGroup{2,2}, grid :: TwoDGrid, 
                                     spline_degree :: Int, smoothing_type  )

        npart = pg.n_particles
        spline1 = SplinePP(spline_degree, grid.nx)
        spline2 = SplinePP(spline_degree, grid.ny)

        n_span = spline_degree + 1

        if smoothing_type == :collocation)
           scaling = 1.0/( grid.dx * grid.dy )
        elseif (smoothing_type == :galerkin)
           scaling = 1.0
        else
           println( "Smoothing Type $smoothing_type not implemented for kernel_smoother_spline_2d. ")
        end

        spline_val = zeros( n_span, 2)
    
        new( grid,  npart, [spline1, spline2] )

    end

end

"""
   compute_shape_factor_spline_2d(self, position, indices)

Helper function computing shape factor
- pm : kernel smoother object
"""
function compute_shape_factor(pm :: ParticleMeshCoupling2D, position, indices)

    x1 = position / pm.grid.delta_x
    indices = ceiling.([xi(1:2))
    xi(1:2) = xi(1:2) - real(indices -1,f64)
    indices =  indices - self%spline_degree
    uniform_bsplines_eval_basis(pm.spline_degree, x1, pm.spline_val(1:self%n_span,1))
    uniform_bsplines_eval_basis(pm.spline_degree, x2, pm.spline_val(1:self%n_span,2))

end 

"""
    add_charge_single_spline_pp_2d(self, position, marker_charge, rho_dofs)

## Add charge of single particle

- Information about the 2d mesh
  * delta_x(2)  : Value of grid spacing along both directions.
  *  domain(2,2) : Definition of the domain: domain(1,1) = x1_min, domain(2,1) = x2_min,  domain(1,2) = x1_max, domain(2,2) = x2_max
- Information about the particles
  * no_particles : Number of particles of underlying PIC method (processor local)
  * n_span : Number of intervals where spline non zero (spline_degree + 1)
  * scaling
  
- position : Particle position
- marker_charge : Particle weight times charge
- rho_dofs : spline coefficient of accumulated density
    
"""
function add_charge_single_spline_pp_2d(self, position, marker_charge, rho_dofs)

    xi(1:2) = (position(1:2) - self%domain(:,1)) /self%delta_x
    indices = floor(xi(1:2))+1
    xi(1:2) = xi(1:2) - real(indices -1,f64)
    indices =  indices - self%spline_degree
    
    spline_pp_horner_m_2d(self%spline_pp, self%spline_val,[self%spline_degree,self%spline_degree], xi)

    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          rho_dofs(index2d) = rho_dofs(index2d) +&
               ( marker_charge* self%scaling * &
               self%spline_val(i1,1) * self%spline_val(i2,2))
       end
    end

end 

"""
    add_charge_single_spline_2d(self, position, marker_charge, rho_dofs)

Add charge of single particle
- position : Particle position
- marker_charge : Particle weight times charge
- rho_dofs : spline coefficient of accumulated density
"""
function add_charge_single_spline_2d(self, position, marker_charge, rho_dofs)
    
    call compute_shape_factor_spline_2d(self, position, indices)
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          rho_dofs(index2d) = rho_dofs(index2d) +&
               ( marker_charge* self%scaling * &
               self%spline_val(i1,1) * self%spline_val(i2,2))
       end
    end

end


"""
     add_current_spline_2d( self, position_old, position_new, marker_charge, j_dofs )

## Add current with integration over x

- position_old : Position of the particle
- position_new : Position of the particle
- marker_charge : Particle weights time charge
- j_dofs : Coefficient vector of the current density

"""
function add_current_spline_2d( self, position_old, position_new, marker_charge, j_dofs )

    print("Running function at ",$("$(__source__.file)"),":",$("$(__source__.line)"))

    error("unimplemented")
    
end 
  

"""
    add_current_update_v_spline_2d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)

## Add current and update v for single particle
- position_old(self%dim) !< Position at time t
- position_new(self%dim) !< Position at time t+\Delta t
- marker_charge !< Particle weight time charge
- qoverm !< Charge over mass ratio
- bfield_dofs(self%n_dofs) !< Coefficient of B-field expansion
- vi(:) !< Velocity of the particle
- j_dofs(self%n_dofs) !< Coefficient of current expansion

"""
function add_current_update_v_spline_2d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)

    error("unimplemented")

end 

  
"""
    add_particle_mass_spline_2d(self, position, marker_charge, particle_mass) 
"""
function add_particle_mass_spline_2d(self, position, marker_charge, particle_mass) 
    class (sll_t_particle_mesh_coupling_spline_2d), intent( inout ) :: self !< Kernel smoother object
    sll_real64,                    intent( in )    :: position(self%dim) !< Position of the particle
    sll_real64,                    intent( in )    :: marker_charge !< Particle weight times charge
    sll_real64,                    intent( inout ) :: particle_mass(:,:) !< Coefficient vector of the charge distribution
    
    error("unimplemented")
    
end 

     
"""
    add_current_update_v_spline_pp_2d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)

Add current and update v for single particle
- position_old(self%dim) !< Position at time t
- position_new(self%dim) !< Position at time t+\Delta t
- marker_charge !< Particle weight time charge
- qoverm !< Charge over mass ratio
- bfield_dofs(self%n_dofs) !< Coefficient of B-field expansion
- vi(:) !< Velocity of the particle
- j_dofs(self%n_dofs) !< Coefficient of current expansion
"""
function add_current_update_v_spline_pp_2d (self, position_old, position_new, marker_charge, qoverm, bfield_dofs, vi, j_dofs)

    error("unimplemented")

end

"""
    evaluate_field_single_spline_pp_2d(self, position, field_dofs_pp, field_value)

Evaluate field at at position \a position using horner scheme

"""
function evaluate_field_single_spline_pp_2d(self, position, field_dofs_pp, field_value)
    class( sll_t_particle_mesh_coupling_spline_2d), intent(inout)  :: self !< kernel smoother object    
    sll_real64,                              intent( in )   :: position(self%dim) !< Position where to evaluate
    sll_real64,                              intent(in)     :: field_dofs_pp(:,:) !< Degrees of freedom in kernel representation.
    sll_real64,                              intent(out)    :: field_value !< Value of the field
       
    !local variables
    sll_real64 :: xi(2)
    sll_int32  :: indices(2)
   
    xi(1:2) = (position(1:2) - self%domain(:,1)) /self%delta_x
    indices = floor(xi(1:2))+1
    xi(1:2) = xi(1:2) - real(indices -1,f64)
     
    field_value = sll_f_spline_pp_horner_2d([self%spline_degree,self%spline_degree], field_dofs_pp, xi, indices,self%n_grid)
        
end subroutine evaluate_field_single_spline_pp_2d
   

"""
    evaluate_field_single_spline_2d(self, position, field_dofs, field_value)

- position(self%dim) : Position where to evaluate
- field_dofs(self%n_dofs) : Degrees of freedom in kernel representation.
- field_value : Value of the field

Evaluate field with given dofs at position \a position
"""
function evaluate_field_single_spline_2d(self, position, field_dofs, field_value)
    
    compute_shape_factor_spline_2d(self, position, indices)

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          field_value = field_value + &
               field_dofs(index2d) *  &
               self%spline_val(i1,1) *&
               self%spline_val(i2,2)
       end
    end

end

"""
    evaluate_multiple_spline_2d(self, position, components, field_dofs, field_value)

## Evaluate multiple fields at position \a position
- position(self%dim) : Position where to evaluate
- components(:) : Components of the field that shall be evaluated
- field_dofs(:,:) : Degrees of freedom in kernel representation.
- field_value(:) : Value of the field
"""
function evaluate_multiple_spline_2d(self, position, components, field_dofs, field_value)
    
    compute_shape_factor_spline_2d(self, position, indices)

    field_value = 0.0_f64
    do i1 = 1, self%n_span
       index1d(1) = indices(1)+i1-2
       do i2 = 1, self%n_span
          index1d(2) = indices(2)+i2-2
          index2d = index_1dto2d_column_major(self,index1d)
          field_value = field_value + &
               field_dofs(index2d,components) *  &
               self%spline_val(i1,1) *&
               self%spline_val(i2,2)
       end
    end

end


"""
    index_1dto2d_column_major(self, index1d) 

Self function computes the index of a 1D array that stores 2D data in column major ordering. It also takes periodic boundary conditions into account.
- index1d(2) !< 2d array with indices along each of the two directions (start counting with zero).
- index2d    !< Corresponding index in 1d array representing 2d data (start counting with one).
"""
function index_1dto2d_column_major(self, index1d)

    index1d(1) = modulo(index1d(1), self%n_grid(1))
    index1d(2) = modulo(index1d(2), self%n_grid(2))
    index2d = index1d(1) + index1d(2)*self%n_grid(1) + 1

    return index2d

end 
