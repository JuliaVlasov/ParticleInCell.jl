using Plots

"""
    vm2d(nstep)

"""
function run( nstep :: Int, time :: Vector{Float64}, energy :: Vector{Float64} )

    nstep = Int32(nstep)

    ccall((:vm2d2v, "./libpic.dylib"), Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}), nstep, time, energy)

    plot( time, energy)

end

nstep = 250
time = zeros(nstep+1)
energy = zeros(nstep+1)

@time run( nstep, time, energy )
