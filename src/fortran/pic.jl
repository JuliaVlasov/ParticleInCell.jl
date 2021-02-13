"""
    vm2d(nstep)

"""
function vm2d( nstep )

    nstep = Int32(nstep) # could overflows

    ccall((:run, "./libpic.dylib"), Cvoid, (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}), size_x, size_y, src, dest, err)
    return err[]

end

using Test
nstep = 250

err = heat(a, b)

@test b == [ 1 1 1; 1 -4 1 ; 1 1 1]
@test err ==  16.0 # ! equality between reals

