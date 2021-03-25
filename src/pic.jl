
module F90

using Libdl
import ..ParticleInCell: TwoDGrid, CloudInCell

const piclib = joinpath(@__DIR__, "libpic")

open(joinpath(@__DIR__, "pic.f90")) do f90file
    open(
        `gfortran -fPIC -w -O3 -shared -x f95 -o $(piclib * "." * Libdl.dlext) -`,
        "w",
    ) do f
        print(f, read(f90file, String))
    end
end



function compute_current!(jx, jy, m::TwoDGrid, kernel :: CloudInCell, p )

    nx = Int32(m.nx)
    ny = Int32(m.ny)
    dx = m.dx
    dy = m.dy
    nbpart = Int32(size(p.array)[2])

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
        p.array,
        jx,
        jy,
    )

    for i = 1:nx+1
        jx[i, 1] += jx[i, ny+1]
        jx[i, ny+1] = jx[i, 1]
    end
    for j = 1:ny+1
        jx[1, j] += jx[nx+1, j]
        jx[nx+1, j] = jx[1, j]
    end

end

function push_x!(p, mesh::TwoDGrid, dt::Float64)

    nbpart = Int32(size(p.array)[2])
    dimx = mesh.dimx
    dimy = mesh.dimy

    ccall(
        (:push_x, piclib),
        Cvoid,
        (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}),
        nbpart,
        dimx,
        dimy,
        p.array,
        dt,
    )

end

export push_v!

function push_v!(p, kernel::CloudInCell, m::TwoDGrid, ex, ey, bz, dt::Float64)

    nbpart = Int32(size(p.array)[2])
    nx = Int32(m.nx)
    ny = Int32(m.ny)
    dx = m.dx
    dy = m.dy

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64},
            Ref{Int32},
            Ref{Int32},
            Ref{Float64},
            Ref{Float64},
            Ptr{Float64},
            Ptr{Float64},
            Ptr{Float64}), nbpart, p.array, dt, nx, ny, dx, dy, ex, ey, bz)

end

end
