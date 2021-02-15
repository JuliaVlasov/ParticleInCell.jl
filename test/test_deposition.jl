using Random

@testset "deposition" begin

dimx, dimy = 4π, 4π
nx, ny = 128, 128
mesh = Mesh( dimx, nx, dimy, ny)
fdtd = FDTD( mesh )
rng = MersenneTwister(1234)
nbpart = 1_000_000
p = zeros(7,nbpart)
randn!(rng, p)
p[1:2,:] .+= 2π
p[3:4,:] .= (dimx * dimy) / nbpart

@test sum(view(p,1,:)) / nbpart ≈ 2π atol=1e-2
@test sum(view(p,2,:)) / nbpart ≈ 2π atol=1e-2

compute_current!( fdtd, p, nbpart)

@test sum(view(fdtd.ebj,4,:,:)) ≈ 0.02493672730470462
@test sum(view(fdtd.ebj,5,:,:)) ≈ 0.02493672730470462

f90_deposition!( fdtd, p)

@test sum(view(fdtd.ebj,4,:,:)) ≈ 0.02493672730470462
@test sum(view(fdtd.ebj,5,:,:)) ≈ 0.02493672730470462

end
