"""
    heat(src::Array{Float64,2}, dest::Array{Float64,2})

apply a stencil to evaluate laplacian on a cartesian grid
```
dest[i,j] = 4 * src[i,j] - src[i-1, j] - src[i+1, j] - src[i, j-1] - src[i, j+1]
```
"""
function heat(src::Array{Float64,2}, dest::Array{Float64,2})
    @assert size(dest) == size(src)
    (size_x, size_y) = size(dest)
    size_x = Int32(size_x) # could overflows
    size_y = Int32(size_y)
    err = Ref{Float64}(0.)
    ccall((:heatKernel, "./libheatKernel.dylib"), Cvoid, (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}), size_x, size_y, src, dest, err)
    return err[]
end

using Test
a = ones(3,3) ; a[2,2] = 0.; b = copy(a)
err = heat(a, b)
@test b == [ 1 1 1; 1 -4 1 ; 1 1 1]
@test err ==  16.0 # ! equality between reals

