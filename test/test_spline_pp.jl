@testset " spline pp 1d " begin

import ParticleInCell: SplinePP, b_to_pp, uniform_bsplines_eval_basis, horner_1d

ncells = 8
  
domain    = [0., 2π]
delta_x   = (domain[2] - domain[1]) / ncells
 
b_coeffs  = rand(ncells)

for degree = 1:5

    xp = 4.8141437173987800
     
    spline_pp = SplinePP( degree, ncells)
    
    pp_coeffs = b_to_pp(spline_pp, ncells, b_coeffs)
    
    xi    = (xp - domain[1])/delta_x
    index = floor(Int64, xi)+1
    xi    = xi - (index-1)
    res   = horner_1d(degree, pp_coeffs, xi, index)
        
    index = index - degree
        
    val = uniform_bsplines_eval_basis(degree, xi)
    
    res2 = 0.0
    
    for i = 1:degree+1 
       index1d = (index+i-2) % ncells + 1
       res2    += b_coeffs[index1d] * val[i]
    end 
    
    @test abs(res-res2) < 1e-15
    
    @test spline_pp.degree - degree < 1e-15
    
    xp    = rand()
    index = 1
    res = horner_1d(degree, pp_coeffs, xp, index)
    
    res2 = 0.
    for i=1:degree+1
        res2 += pp_coeffs[i,1] * xp^((degree+1)-i)
    end
    
    @test abs(res-res2) < 1e-12

end
     
end


@testset " spline pp 2d " begin

    using Random

    ncells = 50
    degree = 3

    b_coeffs = zeros(ncells*ncells)
    pp_coeffs = zeros((degree+1)*(degree+1),ncells*ncells)
       
    dx = 2π / ncells
    dy = 2π / ncells

    rng = MersenneTwister(42)
    rand!(rng, b_coeffs)
    
    spline1 = SplinePP( degree, ncells)
    spline2 = SplinePP( degree, ncells)

    ParticleInCell.b_to_pp_2d!(pp_coeffs, spline1, spline2, b_coeffs)

    xp, yp = rand(rng, 2) .* 2π

    xi = xp / dx
    yi = yp / dy

    ind_x = floor(Int,xi)+1
    ind_y = floor(Int,yi)+1

    xi -= (ind_x-1)
    yi -= (ind_y-1)

    res1 = ParticleInCell.horner_2d((degree, degree), pp_coeffs, 
          (xi, yi), (ind_x, ind_y), (ncells, ncells))

    @show res1
    
    ind_x -= degree
    ind_y -= degree
    
    val1 = uniform_bsplines_eval_basis(degree, xi)
    val2 = uniform_bsplines_eval_basis(degree, yi)
    
    res2 = 0.0
    for i = 1:degree+1
       idx1 = mod(ind_x+i-2, ncells) 
       for j = 1:degree+1
          idx2 = mod( ind_y+j-2, ncells)
          index2d = idx1 + idx2*ncells + 1
          res2 += b_coeffs[index2d] * val1[i] * val2[j]
       end
    end

    @show res2
    @test abs(res1-res2) < 1e-15
    
    #call random_seed()
    #call random_number(xp)
    #res=sll_f_spline_pp_horner_2d(degree, pp_coeffs, xp,[1,1],[1,1])
    #res2=0._f64
    #do i=1, degree(1)+1
    #   do j=1, degree(2)+1
    #      res2=res2+pp_coeffs(i+(j-1)*(degree(1)+1),1)*xp(1)**((degree(1)+1)-i)*xp(2)**((degree(2)+1)-j)
    #   end do
    #end do
    #if(abs(res-res2)>1E-12) then
    #   fail=.true.
    #   print*, xp
    #   print*,'error in horner'
    #end if

    #!print*, 'res-res2=', res-res2
  
end
     
  
