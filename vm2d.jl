using ProgressMeter
using TimerOutputs
using ParticleInCell

include("run_function.jl")

t, energy = run( 1 ) # trigger building
@show nstep = 1000
@time t, energy = run( nstep, npm = 500 )
show(to)

open("results_jl.dat", "w") do f

    for i in 1:nstep
        println(f, t[i], " ", energy[i])
    end

end

println()
