
module F90

using Libdl
import ..ParticleInCell: TwoDGrid

const piclib = joinpath(@__DIR__, "libpic")

open(joinpath(@__DIR__, "pic.f90")) do f90file
    open(
        `gfortran -fPIC -w -O3 -shared -x f95 -o $(piclib * "." * Libdl.dlext) -`,
        "w",
    ) do f
        print(f, read(f90file, String))
    end
end

function interpolation!(p::Array{Float64,2}, m::TwoDGrid)

    nx = Int32(m.nx)
    ny = Int32(m.ny)
    dx = m.dx
    dy = m.dy
    nbpart = Int32(size(p)[2])

    ccall(
        (:interpolation, piclib),
        Cvoid,
        (
            Ref{Int32},
            Ref{Int32},
            Ref{Int32},
            Ref{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
        ),
        nbpart,
        nx,
        ny,
        dx,
        dy,
        m.ex,
        m.ey,
        m.bz,
        p,
    )

end

function compute_current!(m::TwoDGrid, p::Array{Float64,2})

    nx = Int32(m.nx)
    ny = Int32(m.ny)
    dx = m.dx
    dy = m.dy
    nbpart = Int32(size(p)[2])

    ccall(
        (:deposition, piclib),
        Cvoid,
        (
            Ref{Int32},
            Ref{Int32},
            Ref{Int32},
            Ref{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64},
        ),
        nbpart,
        nx,
        ny,
        dx,
        dy,
        p,
        m.jx,
        m.jy,
    )

    for i = 1:nx+1
        m.jx[i, 1] += m.jx[i, ny+1]
        m.jx[i, ny+1] = m.jx[i, 1]
    end
    for j = 1:ny+1
        m.jx[1, j] += m.jx[nx+1, j]
        m.jx[nx+1, j] = m.jx[1, j]
    end

end

function push_x!(p::Array{Float64,2}, mesh::TwoDGrid, dt::Float64)

    nbpart = Int32(size(p)[2])
    dimx = mesh.dimx
    dimy = mesh.dimy

    ccall(
        (:push_x, piclib),
        Cvoid,
        (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}),
        nbpart,
        dimx,
        dimy,
        p,
        dt,
    )

end

export push_v!

function push_v!(p::Array{Float64,2}, dt::Float64)

    nbpart = Int32(size(p)[2])

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64}), nbpart, p, dt)

end

end
