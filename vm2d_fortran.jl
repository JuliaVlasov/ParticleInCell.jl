using LinearAlgebra, ProgressMeter
using TimerOutputs, Sobol

import ParticleInCell: Mesh, landau_sampling!, FDTD, faraday!, ampere_maxwell!
import ParticleInCell: compute_energy, update_fields!
import ParticleInCell.F90: interpolation!, compute_current!, push_x!, push_v!

include("run_function.jl")

t, energy = run( 1 )
@show nstep = 1000
@time t, energy = run( nstep, npm=500 )
show(to)

open("results_f90.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
