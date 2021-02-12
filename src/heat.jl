function heat(src::Array{Float64,2}, dest::Array{Float64,2})
    @assert size(dest) == size(src)
    (size_x, size_y) = size(dest)
    size_x = Int32(size_x)
    size_y = Int32(size_y)
    err = Ref{Float64}(0.) 
    ccall((:heatKernel, "./libheatKernel.so"), Cvoid, (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ref{Float64}), size_x, size_y, src, dest, err)
    return err[]
end
