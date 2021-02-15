using Random

@testset "interpolation" begin

dimx, dimy = 4π, 4π
nx, ny = 128, 128
mesh = Mesh( dimx, nx, dimy, ny)
fdtd = FDTD( mesh )
rng = MersenneTwister(1234)
nbpart = 1_000_000
p = zeros(7,nbpart)
randn!(rng, p)
p[1:2,:] .+= 2π
p[3:7,:] .= 0

for i in 1:nx+1, j in 1:ny+1
    x =  mesh.x[i] - 2π
    y =  mesh.y[j] - 2π
    fdtd.ebj[1,i,j] = exp.(-0.5 * ( x^2 + y^2 ))
    fdtd.ebj[2,i,j] = exp.(-0.5 * ( x^2 + y^2 ))
    fdtd.ebj[3,i,j] = exp.(-0.5 * ( x^2 + y^2 ))
end

interpol_eb!( p, nbpart, fdtd)

@test sum(view(p,5,:)) ≈ 499867.32386298594
@test sum(view(p,6,:)) ≈ 499867.32386298594 
@test sum(view(p,7,:)) ≈ 499867.32386298594 

f90_interpolation!( p, fdtd)

@test sum(view(p,5,:)) ≈ 499867.32386298594 
@test sum(view(p,6,:)) ≈ 499867.32386298594 
@test sum(view(p,7,:)) ≈ 499867.32386298594 

end
