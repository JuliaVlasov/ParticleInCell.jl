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
    spline1 :: SplinePP
    spline2 :: SplinePP
    n_span :: Int
    spline_degree :: Int
    scaling :: Float64
    spline_val :: Array{Float64,2}

    function ParticleMeshCoupling2D( pg :: ParticleGroup{2,2}, grid :: TwoDGrid, 
                                     spline_degree :: Int, smoothing_type  )

        npart = pg.n_particles
        spline1 = SplinePP(spline_degree, grid.nx)
        spline2 = SplinePP(spline_degree, grid.ny)

        n_span = spline_degree + 1

        if smoothing_type == :collocation
           scaling = 1.0/( grid.dx * grid.dy )
        elseif smoothing_type == :galerkin
           scaling = 1.0
        else
           println( "Smoothing Type $smoothing_type not implemented for kernel_smoother_spline_2d. ")
        end

        spline_val = zeros( n_span, 2)
    
        new( grid, npart, spline1, spline2, n_span, spline_degree, scaling, spline_val )

    end

end


"""
   compute_shape_factor_spline_2d(self, position, indices)

Helper function computing shape factor
- pm : kernel smoother object
"""
function compute_shape_factor!(pm :: ParticleMeshCoupling2D, xp, yp)

    xp = (xp - pm.grid.xmin) / pm.grid.dx
    yp = (yp - pm.grid.ymin) / pm.grid.dy
    ip = ceil(Int, xp)
    jp = ceil(Int, yp)
    dxp = xp - (ip-1)
    dyp = yp - (jp-1)

    uniform_bsplines_eval_basis!( pm.spline_val, pm.spline_degree, dxp, dyp) 

    return (ip - pm.spline_degree, jp - pm.spline_degree)

end 

"""
    index_1dto2d_column_major(pm, index1d) 

Self function computes the index of a 1D array that stores 2D data in column major ordering. 
It also takes periodic boundary conditions into account.
- index1d_1 !< indice along x (start counting with zero).
- index1d_2 !< indice along y (start counting with zero).
- index2d   !< Corresponding index in 1d array representing 2d data (start counting with one).
"""
function index_1dto2d_column_major(pm, index1d_1, index1d_2)

    index1d_1 = mod(index1d_1, pm.grid.nx)
    index1d_2 = mod(index1d_2, pm.grid.ny)
    index2d = index1d_1 + index1d_2 * pm.grid.nx + 1

    return index2d

end 


export add_charge!

"""
    add_charge(self, x_position, y_position marker_charge, rho_dofs)

Add charge of single particle
- position : Particle position
- marker_charge : Particle weight times charge
- rho_dofs : spline coefficient of accumulated density
"""
function add_charge!(rho_dofs, pm ::  ParticleMeshCoupling2D , xp, yp, marker_charge)
    
    ind_x, ind_y = compute_shape_factor!(pm, xp, yp )

    for i1 = 1:pm.n_span
       index1d_1 = ind_x + i1 - 2
       for i2 = 1:pm.n_span
          index1d_2 = ind_y + i2 -2
          index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
          @show pm.spline_val[i1,1], pm.spline_val[i2,2], marker_charge, pm.scaling
          rho_dofs[index2d] += ( marker_charge * pm.scaling * pm.spline_val[i1,1] * pm.spline_val[i2,2])
       end
    end

end


export add_charge_pp!
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
function add_charge_pp!(rho_dofs, pm :: ParticleMeshCoupling2D, xp, yp, marker_charge)

    xp = (xp - pm.grid.xmin) / pm.grid.dx
    yp = (yp - pm.grid.ymin) / pm.grid.dy
    ip = floor(Int, xp)+1
    jp = floor(Int, yp)+1
    dxp = xp - (ip-1)
    dyp = yp - (jp-1)

    ip =  ip - pm.spline_degree
    jp =  jp - pm.spline_degree
    
    horner_m_2d!(pm.spline1, pm.spline2, pm.spline_val, pm.spline_degree, dxp, dyp)

    for i1 = 1:pm.n_span
       index1d_1 = ip+i1-2
       for i2 = 1:pm.n_span
          index1d_2 = jp+i2-2
          index2d = index_1dto2d_column_major(pm, index1d_1, index1d_2)
          rho_dofs[index2d] += ( marker_charge * pm.scaling * pm.spline_val[i1,1] * pm.spline_val[i2,2])
       end
    end

end 


#=

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


=#
