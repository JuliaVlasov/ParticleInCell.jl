using Libdl

const piclib = joinpath(@__DIR__, "libpic.dylib")

open(joinpath(@__DIR__, "pic.f90")) do f90file
    open(`gfortran -fPIC -w -O3 -shared -x f95 -o $(piclib * "." * Libdl.dlext) -`, "w") do f
        print(f, read(f90file, String))
    end
end

export vm2d2v

function vm2d2v( nstep :: Int, time :: Vector{Float64}, energy :: Vector{Float64} )

    nstep = Int32(nstep)

    ccall((:vm2d2v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}), nstep, time, energy)

end

export interpolation!

function interpolation!( p, fdtd :: FDTD )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(particles)[2])

    f = fdtd.ebj

    ccall((:interpolation, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, f, p)

end

export deposition!

function deposition!( fdtd :: FDTD, p )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(particles)[2])

    f = fdtd.ebj
    jx = fdtd.jx
    jy = fdtd.jy

    ccall((:deposition, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, p, f, jx, jy)

end

function f90_push_x!( p, nbpart, dimx, dimy, dt)

    nbpart = Int32(nbpart)

    ccall((:push_x, piclib), Cvoid, (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}), nbpart, dimx, dimy, p.data, dt)

end

function f90_push_v!( p, nbpart, dt :: Float64)

    nbpart = Int32(nbpart)

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64}), nbpart, p, dt)

end
