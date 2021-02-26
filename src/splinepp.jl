using StaticArrays

const  inv_2   = 1. / 2.
const  inv_3   = 1. / 3.
const  inv_4   = 1. / 4.
const  inv_6   = 1. / 6.
const  inv_8   = 1. / 8.
const  inv_12  = 1. / 12.
const  inv_18  = 1. / 18.
const  inv_20  = 1. / 20.
const  inv_24  = 1. / 24.
const  inv_30  = 1. / 30.
const  inv_36  = 1. / 36.
const  inv_48  = 1. / 48.
const  inv_72  = 1. / 72.
const  inv_120 = 1. / 120.
const  inv_144 = 1. / 144.
const  inv_720 = 1. / 720.

"""
    SplinePP( degree, ncells)

- `degree` : degree of 1d spline
- `poly_coeffs` : `poly_coeffs[i,j]` coefficient of ``x^{deg+1-j}`` for ith B-spline function  size= (degree+1, degree+1)
- `poly_coeffs_fp` : `poly_coeffs[i,j]` coefficient of ``x^{deg+1-j}`` for ith B-spline function  size= (degree+1, degree+1)
- `ncells` : number of gridcells
- `scratch_b` : scratch data for `b_to_pp-converting`
- `scratch_p` : scratch data for `b_to_pp-converting`
"""
struct SplinePP{N}

    degree         :: Int
    poly_coeffs    :: SMatrix{N,N,Float64}
    poly_coeffs_fp :: SMatrix{N,N,Float64}
    ncells         :: Int

    function SplinePP(degree, ncells) 
    
        @assert (ncells >= degree )
    
        if degree == 1

             poly_coeffs    = SMatrix{2,2}( -1., 1., 1., 0. ) 
             poly_coeffs_fp = SMatrix{2,2}( -inv_2, 1., inv_2, 0. )
           
        elseif degree == 2

             poly_coeffs    = SMatrix{3,3}(inv_2, -1., inv_2 , -1., 1., 
                                       inv_2, inv_2, 0., 0.)
             poly_coeffs_fp = SMatrix{3,3}(inv_6, -inv_2, inv_2 , -inv_3, 
                                       inv_2, inv_2, inv_6, 0., 0)

        elseif degree == 3

            poly_coeffs = SMatrix{4,4}( -inv_6, inv_2, -inv_2, inv_6
                                  ,  inv_2,   -1.,     0., 4*inv_6
                                  , -inv_2, inv_2,  inv_2, inv_6
                                  ,  inv_6,    0.,     0., 0.)
           
            poly_coeffs_fp = SMatrix{4,4}(- inv_24, inv_6, -inv_4, inv_6
                                     ,  inv_8, -inv_3,     0., 4*inv_6
                                     , -inv_8,  inv_6,  inv_4, inv_6
                                     ,  inv_24,    0.,     0., 0.)

        elseif degree == 4

            poly_coeffs = SMatrix{5,5}(inv_24,-inv_6, inv_4,-inv_6, inv_24
                                 ,- inv_6, inv_2,-inv_4,-inv_2, 11*inv_24
                                 ,  inv_4,-inv_2,-inv_4, inv_2, 11*inv_24
                                 ,- inv_6, inv_6, inv_4, inv_6, inv_24
                                 , inv_24,    0.,    0.,    0., 0.   )

            poly_coeffs_fp = SMatrix{5,5}( inv_120,- inv_24, inv_12,-inv_12, inv_24
                                    , - inv_30,  inv_8,-inv_12,-inv_2, 11*inv_24
                                    ,   inv_20,- inv_8,-inv_12,inv_4,11*inv_24
                                    , - inv_30,  inv_24,inv_12,inv_12,inv_24
                                    ,   inv_120,     0.,    0.,    0., 0.)

        elseif degree == 5

            poly_coeffs = SMatrix{6,6}(-inv_120,inv_24,-inv_12,inv_12,-inv_24,inv_120
                                   ,inv_24,-inv_6,inv_6,inv_6,-5*inv_12, 26*inv_120 
                                   ,-inv_12,inv_4,0.,-inv_2,0.,11*inv_20
                                   ,inv_12,-inv_6,-inv_6,inv_6,5*inv_12,26*inv_120 
                                   ,-inv_24,inv_24,inv_12,inv_12,inv_24,inv_120
                                   ,inv_120,0.,0.,0.,0.,0.)

            poly_coeffs_fp = SMatrix{6,6}(-inv_720,inv_120,-inv_48,inv_36,-inv_48,inv_120
                                      , inv_144,-inv_30,inv_24,inv_18,-5*inv_24, 26*inv_120
                                      ,-inv_72,inv_20,0.,-inv_6,0.,11*inv_20
                                      ,inv_72,-inv_30,-inv_24,inv_18,5*inv_24,26*inv_120
                                      ,-inv_144,inv_120,inv_48,inv_36,inv_48,inv_120
                                      ,inv_720,0.,0.,0.,0.,0.) 
        else

           throw(ArgumentError(" degree $degree not implemented"))

        end

        N = degree+1

        new{N}( degree, poly_coeffs, poly_coeffs_fp, ncells)

    end 
     
end 

"""
    b_to_pp( SplinePP, ncells, b_coeffs)

Convert 1d spline in B form to spline in pp form with 
periodic boundary conditions
"""
function b_to_pp( self :: SplinePP, ncells :: Int, b_coeffs :: Vector{Float64})

    degp1     = self.degree+1
    pp_coeffs = zeros(Float64, (degp1,ncells)) 
    coeffs    = zeros(Float64, degp1) 
       
    @inbounds for i=1:self.degree   
        coeffs .= vcat(b_coeffs[end-self.degree+i:end],b_coeffs[1:i]) 
        for j=1:degp1
            pp_coeffs[j, i] = sum(coeffs .* self.poly_coeffs[j,:])
        end
    end
    
    @inbounds for i=self.degree+1:ncells
        coeffs .= b_coeffs[i-self.degree:i]
        for j=1:degp1
            pp_coeffs[j, i] = sum(coeffs .* self.poly_coeffs[j,:])
        end
    end

    pp_coeffs

end

"""
    horner_1d(degree, pp_coeffs, x, index)

Perform a 1d Horner schema on the `pp_coeffs` at index
"""
function horner_1d(degree :: Int, pp_coeffs, x :: Float64, index :: Int)
    
    res = pp_coeffs[1,index]
    @inbounds for i=1:degree
       res = res * x + pp_coeffs[i+1,index]
    end
    res

end

"""
    horner_primitive_1d(val, degree, pp_coeffs, x)

Perform a 1d Horner schema on the `pp_coeffs` evaluate at x
"""
function horner_primitive_1d(val :: Vector{Float64}, degree, pp_coeffs, x)

  @inbounds for i in eachindex(val)

     val[i] = horner_1d(degree, pp_coeffs, x, i) * x

  end

end 


"""
    b_to_pp_2d_periodic( spl, n_cells, b_coeffs, pp_coeffs)

Convert 2d spline in B form to spline in pp form with periodic boundary 
conditions. This is a special case of the procedure spline_pp_b_to_pp_2d 
for the double periodic case to avoid the select case statements
- spl : arbitrary degree 2d spline 
- n_cells : number of gridcells
- b_coeffs(n_cells(1)*n_cells(2)) : coefficients of spline in B-form
- pp_coeffs((spl%spline1%degree+1)*(spl%spline2%degree+1),n_cells(1)*n_cells(2)) : coefficients of spline in pp-form
"""
function b_to_pp_2d_periodic( spl, n_cells, b_coeffs, pp_coeffs)

    degree1 = spl[1].degree
    degree2 = spl[2].degree

    for j=1:n_cells(2), i=1:n_cells(1)
        b_to_pp_2d_cell(spl, n_cells, b_coeffs, pp_coeffs, i,j)
    end

end 

"""
    horner_m_2d(spl, val, degree, x)

Perform two times a 1d hornerschema on the poly_coeffs
- val : array of values
- degree : degree of the spline
- x : point at which we evaluate our spline
"""
function horner_m_2d!(spline1, spline2, val, degree, x, y)

    for i=1:degree+1
       val[i,1] = horner_1d(degree, spline1.poly_coeffs, x, i)
       val[i,2] = horner_1d(degree, spline2.poly_coeffs, y, i)
    end

end 

#=
"""
    spline_pp_b_to_pp_2d( spl, n_cells, b_coeffs, pp_coeffs)

Convert 2d spline in B form to spline in pp form   
- n_cells(2) : number of gridcells
- b_coeffs   : coefficients of spline in B-form
- pp_coeffs  : coefficients of spline in pp-form
"""
function spline_pp_b_to_pp_2d( spl, n_cells, b_coeffs, pp_coeffs)


    degree1 = spl.spline1.degree
    degree2 = spl.spline2.degree
    n_coeffs1 = spl.spline1.n_coeffs
	n_coeffs2 = spl.spline2.n_coeffs

    select case ( spl%spline1%boundary_conditions )
    case (sll_p_boundary_periodic )
       offset1 = -degree1
    case (sll_p_boundary_clampeddiri:sll_p_boundary_clampeddiri_clamped)
       offset1 = -1
    case default
       offset1 = 0
    end select

    upper(1) = n_cells(1)-degree1+1
    upper(2) = n_cells(2)-degree2+1
    for j=1, n_cells(2)
       for i=1, n_cells(1)
          cell = (j-1)*n_cells(1)+i
          for l2=0,degree2
             select case ( spl%spline2%boundary_conditions )
             case( sll_p_boundary_periodic )
                offset2 = modulo(j-degree2+l2-1,n_coeffs(2))*n_coeffs(1) ! periodic
             case( sll_p_boundary_clampeddiri )
                offset2 = (j+l2-2)*n_coeffs(1)
                if ( j==1 .and. l2 == 0 ) 
                   pp_coeffs_local(:,1) = 0.0_f64
                   cycle
                elseif ( j == n_cells(2) .and. l2 == degree2 ) 
                   pp_coeffs_local(:,degree2+1) = 0.0_f64
                   cycle
                end
             case( sll_p_boundary_clamped_clampeddiri )
                offset2 = (j+l2-1)*n_coeffs(1)
                if ( j == n_cells(2) .and. l2 == degree2 ) 
                   pp_coeffs_local(:,degree2+1) = 0.0_f64
                   cycle
                end
             case( sll_p_boundary_clampeddiri_clamped )
                offset2 = (j+l2-2)*n_coeffs(1)
                if ( j==1 .and. l2 == 0 ) 
                   pp_coeffs_local(:,1) = 0.0_f64
                   cycle
                end
             case default
                offset2 = (j+l2-1)*n_coeffs(1)
             end select

             for l1=0,degree1
                select case( spl.spline1.boundary_conditions )
                case ( sll_p_boundary_periodic )
                   index = modulo(i+offset1+l1-1,n_coeffs(1))+1
                   !print*, i, j, l1, l2, offset2, index
                   spl.spline1.scratch_b(l1+1) = b_coeffs( offset2+index )
                case ( sll_p_boundary_clamped )
                   index = i + offset1 + l1
                   spl.spline1.scratch_b(l1+1) = b_coeffs( offset2+index )
                   ! Set to zero for Dirichlet if we are at the boundary
                case ( sll_p_boundary_clamped_clampeddiri )
                   if ( i == n_cells(1) .and. l1 == degree1 ) 
                      spl.spline1.scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      spl.spline1.scratch_b(l1+1) = b_coeffs( offset2+index )
                   end
                case ( sll_p_boundary_clampeddiri_clamped )
                   if ( i == 1 .and. l1 == 0 ) 
                      spl.spline1.scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      spl.spline1.scratch_b(l1+1) = b_coeffs( offset2+index )
                   end
                case ( sll_p_boundary_clampeddiri )
                   if ( i == 1 .and. l1 == 0 ) 
                      spl.spline1.scratch_b(l1+1) = 0.0_f64
                   elseif ( i == n_cells(1) .and. l1 == degree1 ) 
                      spl.spline1.scratch_b(l1+1) = 0.0_f64
                   else
                      index = i + offset1 + l1
                      spl.spline1.scratch_b(l1+1) = b_coeffs( offset2+index )
                   end
                end select
             end
             ! For clamped splines, we need to use the boundary coefficients
             select case( spl.spline1.boundary_conditions )
             case ( sll_p_boundary_periodic )
                call sll_s_spline_pp_b_to_pp_1d_cella(degree1, spl.spline1.poly_coeffs, &
                     spl%spline1%scratch_b,pp_coeffs_local(:,l2+1))
             case default
                if ( i > degree1-1 ) 
                   if ( i < upper(1)+1 ) 
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree1, spl%spline1%poly_coeffs, &
                           spl%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                   else
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                           spl%spline1%poly_coeffs_boundary_right(:,:,i-upper(1)), &
                           spl%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                   end
                else

                   call sll_s_spline_pp_b_to_pp_1d_cella(degree1, &
                        spl%spline1%poly_coeffs_boundary_left(:,:,i), &
                        spl%spline1%scratch_b,pp_coeffs_local(:,l2+1))
                end


             end select
          end
          for l1=0,degree1
             spl%spline2%scratch_b =  pp_coeffs_local(l1+1,:)
             select case( spl%spline2%boundary_conditions )
             case ( sll_p_boundary_periodic )
                call sll_s_spline_pp_b_to_pp_1d_cella(degree2, spl%spline2%poly_coeffs, &
                     spl%spline2%scratch_b,spl%spline2%scratch_p)
             case default
                if ( j > degree2-1 ) 
                   if ( j < upper(2)+1 ) 
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                           spl%spline2%poly_coeffs, &
                           spl%spline2%scratch_b,spl%spline2%scratch_p)
                   else
                      call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                           spl%spline2%poly_coeffs_boundary_right(:,:,j-upper(2)), &
                           spl%spline2%scratch_b,spl%spline2%scratch_p)
                   end
                else

                   call sll_s_spline_pp_b_to_pp_1d_cella(degree2, &
                        spl%spline2%poly_coeffs_boundary_left(:,:,j), &
                        spl%spline2%scratch_b,spl%spline2%scratch_p)
                end


             end select


             !call sll_s_spline_pp_b_to_pp_1d_cell(spl%spline2, &
             !     spl%spline2%scratch_b,spl%spline2%scratch_p)
             for l2=0,degree2
                pp_coeffs( l2*(degree1+1)+l1+1 , cell) = spl%spline2%scratch_p(l2+1)
             end
          end
       end
    end

end 


"""
    spline_pp_b_to_pp_2d_cell(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)

Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions 

- spline1 : arbitrary degree 1d spline
- spline2 : arbitrary degree 1d spline
- n_cells(2) : number of gridcells
- b_coeffs(n_cells(1)*n_cells(2)) : coefficients of spline in B-form
- pp_coeffs((spline1.degree+1)*(spline2.degree+1),n_cells(1)*n_cells(2)) : coefficients of spline in pp-form
"""
function spline_pp_b_to_pp_2d_cell(spline1,spline2,n_cells, b_coeffs, pp_coeffs,i,j)

    degree1= spline1.degree
    degree2= spline2.degree
    degp1=degree1+1
    degp2=degree2+1
    # convert b-coefficients in pp-coefficients in first dimension
    if (i>degree1)
       if(j>degree2)
          for l=0:degree2
             spline1.scratch_b=b_coeffs(i-degree1+(j-degp2+l)*n_cells(1):i+(j-degp2+l)*n_cells(1))
             spline_pp_b_to_pp_1d_cell(spline1, spline1.scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end
       else 
          # use of modulo for boundary cells in second dimension 
          for l=0,degree2
             spline1.scratch_b=b_coeffs(i-degree1+modulo(j-degp2+l,n_cells(2))*n_cells(1):i+modulo(j-degp2+l,n_cells(2))*n_cells(1))
             spline_pp_b_to_pp_1d_cell(spline1, spline1.scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
          end
       end
    else 
       !> use of modulo for boundary cells in both dimensions 
       for l=0,degree2
          for k=0, degree1
             spline1.scratch_b(k+1)=b_coeffs(modulo(i-degp1+k,n_cells(1))+1+modulo(j-degp2+l,n_cells(2))*n_cells(1))
          end

          spline_pp_b_to_pp_1d_cell(spline1, spline1.scratch_b, pp_coeffs(1+l*degp1:degp1*(l+1),i+n_cells(1)*(j-1)))
       end
    end

    # convert b-coefficients in pp_coefficients in second dimension
    for k=1:degp1
       for l=1:degp2
          spline2.scratch_b(l)=pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))
       end
       call sll_s_spline_pp_b_to_pp_1d_cell(spline2, spline2.scratch_b, spline2.scratch_p)
       for l=1:degp2
          pp_coeffs(k+degp1*(l-1),i+n_cells(1)*(j-1))=spline2.scratch_p(l)
       end
    end
end 



"""
Perform a 2d hornerschema on the pp_coeffs at the indices
- degree : degree of the spline
- pp_coeffs : coefficients of spline in pp-form
- x : point at which we evaluate our spline
- indices(2) : indices of cell in which is x
- n_cells(2) : number of gridcells
- res : value of the splinefunction at point x
"""
function spline_pp_horner_2d(degree, pp_coeffs, x, indices, n_cells) 
    sll_real64 :: pp_coeffs_1d(degree(2)+1,1)
    sll_int32  :: i
    # Perform a 1d hornerschema in the first dimension
    for i=0:degree(2)
       pp_coeffs_1d(i+1,1)=spline_pp_horner_1d(degree(1), pp_coeffs(1+i*(degree(1)+1):(degree(1)+1)*(i+1),1+(indices(2)-1)*n_cells(1):n_cells(1)*indices(2)), x(1),indices(1))
    end

    #Perform a 1d hornerschema in the second dimension
    spline_pp_horner_1d(degree(2), pp_coeffs_1d, x(2),1)

end 
=#
