using Libdl

const piclib = joinpath(@__DIR__, "libpic")

open(joinpath(@__DIR__, "pic.f90")) do f90file
    open(`gfortran -fPIC -w -O3 -shared -x f95 -o $(piclib * "." * Libdl.dlext) -`, "w") do f
        print(f, read(f90file, String))
    end
end

export f90_interpolation!

function f90_interpolation!( p :: Array{Float64,2}, fdtd :: FDTD )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(p)[2])

    f = fdtd.ebj

    ccall((:interpolation, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, f, p)

end

export f90_deposition!

function f90_deposition!( fdtd :: FDTD, p :: Array{Float64, 2} )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(p)[2])

    f = fdtd.ebj
    jx = fdtd.jx
    jy = fdtd.jy

    ccall((:deposition, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, p, f, jx, jy)

end

function f90_push_x!( p :: Array{Float64,2}, nbpart :: Int, dimx :: Float64, dimy :: Float64, dt :: Float64)

    nbpart = Int32(nbpart)

    ccall((:push_x, piclib), Cvoid, (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}), nbpart, dimx, dimy, p, dt)

end

function f90_push_v!( p :: Array{Float64, 2}, nbpart :: Int, dt :: Float64)

    nbpart = Int32(nbpart)

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64}), nbpart, p, dt)

end
