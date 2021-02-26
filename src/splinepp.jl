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
    horner_m_2d(spl, val, degree, x)

Perform two times a 1d hornerschema on the poly_coeffs
- val : array of values
- degree : degree of the spline
- x : point at which we evaluate our spline
"""
function horner_m_2d!(val, spline1 :: SplinePP, spline2 :: SplinePP, degree, x, y)

    for i=1:degree+1
       val[i,1] = horner_1d(degree, spline1.poly_coeffs, x, i)
       val[i,2] = horner_1d(degree, spline2.poly_coeffs, y, i)
    end

end 



"""
    b_to_pp_2d_cell(spline1, spline2, b_coeffs, pp_coeffs, i, j)

Convert 2d spline in B form in a cell to spline in pp form with periodic boundary conditions 

- spline1 : arbitrary degree 1d spline
- spline2 : arbitrary degree 1d spline
- n_cells(2) : number of gridcells
- b_coeffs(n_cells(1)*n_cells(2)) : coefficients of spline in B-form
- pp_coeffs((spline1.degree+1)*(spline2.degree+1),n_cells(1)*n_cells(2)) : coefficients of spline in pp-form
"""
function b_to_pp_2d_cell!(pp_coeffs, spline1, spline2, b_coeffs, i, j)

    degree1 = spline1.degree
    degree2 = spline2.degree
    degp1 = degree1+1
    degp2 = degree2+1
    ncells1 = spline1.ncells
    ncells2 = spline2.ncells

    # convert b-coefficients in pp-coefficients in first dimension
    if (i>degree1)
        if(j>degree2)
            for l=0:degree2

                scratch_b = view(b_coeffs,i-degree1+(j-degp2+l)*ncells1:i+(j-degp2+l)*ncells1)
                pp_coeffs1 = view(pp_coeffs, 1+l*degp1:degp1*(l+1),i+ncells1*(j-1))

                fill!(pp_coeffs1, 0.0)
                for i=1:degp1, j=1:degp1                                                                        
                    pp_coeffs1[j] += scratch_b[i] * spline1.poly_coeffs[j,i]
                end

            end
        else 
            # use of modulo for boundary cells in second dimension 
            for l=0:degree2
               scratch_b = view(b_coeffs,i-degree1+mod(j-degp2+l,ncells2)*ncells1:i+mod(j-degp2+l,ncells2)*ncells1)
               pp_coeffs1 = view(pp_coeffs,1+l*degp1:degp1*(l+1),i+ncells1*(j-1))

               fill!(pp_coeffs1, 0.0)
               for i=1:degp1, j=1:degp1                                                                        
                   pp_coeffs1[j] += scratch_b[i] * spline.poly_coeffs[j,i]
               end
            end
        end
    else 
        # use of modulo for boundary cells in both dimensions 
        for l=0:degree2
             scratch_b = view(b_coeffs,i-degree1+mod(j-degp2+l,ncells2)*ncells1:i+mod(j-degp2+l,ncells2)*ncells1)
             for k=0:degree1
                scratch_b[k+1]=b_coeffs[mod(i-degp1+k,ncells1)+1+mod(j-degp2+l,ncells2)*ncells1]
             end
             pp_coeffs1 = view(pp_coeffs,1+l*degp1:degp1*(l+1),i+ncells1*(j-1))

             fill!(pp_coeffs, 0.0)
             for i=1:degp1, j=1:degp1                                                                        
                 pp_coeffs[j] += b_coeffs[i] * spline1.poly_coeffs[j,i]
             end
        end
    end

# convert b-coefficients in pp_coefficients in second dimension
    for k=1:degp1
        for l=1:degp2
           spline2.scratch_b[l]=pp_coeffs[k+degp1*(l-1),i+ncells1*(j-1)]
        end
        b_to_pp_1d_cell(spline2, spline2.scratch_b, spline2.scratch_p)

              fill!(pp_coeffs, 0.0)
              for i=1:degp1, j=1:degp1                                                                        
                  pp_coeffs[j] += b_coeffs[i] * spline.poly_coeffs[j,i]
              end

        for l=1:degp2
           pp_coeffs[k+degp1*(l-1),i+ncells1*(j-1)] = spline2.scratch_p[l]
        end

    end
end 

"""
    b_to_pp_2d( spl, n_cells, b_coeffs, pp_coeffs)

Convert 2d spline in B form to spline in pp form   
- n_cells(2) : number of gridcells
- b_coeffs   : coefficients of spline in B-form
- pp_coeffs  : coefficients of spline in pp-form
"""
function b_to_pp_2d!( pp_coeffs, spline1 :: SplinePP, spline2 :: SplinePP, b_coeffs)

    ncells1 = spline1.ncells
    ncells2 = spline2.ncells

    #for j=1:ncells2, i=1:ncells1
    #    b_to_pp_2d_cell!( pp_coeffs, spline1, spline2, b_coeffs, i, j)
    #end

end 

"""
    horner_2d(degree, pp_coeffs, x, indices, ncells) 

Perform a 2d hornerschema on the pp_coeffs at the indices
- degree : degree of the spline
- pp_coeffs : coefficients of spline in pp-form
- x : point at which we evaluate our spline
- indices(2) : indices of cell in which is x
- ncells(2) : number of gridcells
- res : value of the splinefunction at point x
"""
function horner_2d(degrees, pp_coeffs, position, indices, ncells) 

    x1, x2 = position
    d1, d2 = degrees
    i1, i2 = indices
    n1, n2 = ncells

    pp_coeffs_1d = zeros(d2+1,1)

    for i=0:d2
       pp_coeffs_1d[i+1,1] = horner_1d(d1, pp_coeffs[1+i*(d1+1):(d1+1)*(i+1),1+(i2-1)*n1:n1*i2], x1, i1 )
    end

    res = pp_coeffs_1d[1,1]
    for i=1:d2
       res = res * x2 + pp_coeffs_1d[i+1,1]
    end

    res

end 
