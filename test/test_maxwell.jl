using ParticleInCell 
using LinearAlgebra


@testset "Maxwell solver" begin


dimx, dimy = 1, 1
nx, ny = 100, 100
e0 = 1
c = 1
dt = 0.001
nstep = 10

mesh = Mesh( dimx, nx, dimy, ny )

maxwell = MaxwellSolver( mesh, c, e0) 

omega = c * sqrt((2π/dimx)^2+(2π/dimy)^2)

ex = zeros(nx+1,ny+1)
ey = zeros(nx+1,ny+1)
bz = zeros(nx+1,ny+1)

# Ex and Ey are set at t = 0.0
# Bz is set at  t = -dt/2

xn = mesh.x
yn = transpose(mesh.y)

xc = 0.5 .* ( xn[1:end-1] .+ xn[2:end])
yc = transpose(0.5 .* ( yn[1:end-1] .+ yn[2:end]))

@show size(xc)
@show size(yc)

@. maxwell.ex = 0
@. maxwell.ey = 0
@. maxwell.bz = - cos(2π*xc) * cos(2π*yc) * cos(omega*(-0.5*dt))

t = 0

sol_bz = - cos.(2π*xn) .* cos.(2π*yn)
sol_ex = + cos.(2π*xn) .* sin.(2π*yn) .* 2π / omega
sol_ey = - sin.(2π*xn) .* cos.(2π*yn) .* 2π / omega
    
for istep = 1:nstep # Loop over time
    
    faraday!(maxwell, bz, dt)   

    t = t + 0.5dt

    @test maximum(abs.(bz .- sol_bz .* cos(omega*t))) < 0.001

    ampere_maxwell!(maxwell, ex, ey, dt) 

    t = t + 0.5dt

    @test maximum(abs.(ex .- sol_ex .* sin(omega*t))) < 0.001
    @test maximum(abs.(ey .- sol_ey .* sin(omega*t))) < 0.001

    
end # next time step
    
end
