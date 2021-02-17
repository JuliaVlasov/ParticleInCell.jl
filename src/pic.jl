using Libdl

const piclib = joinpath(@__DIR__, "libpic")

open(joinpath(@__DIR__, "pic.f90")) do f90file
    open(`gfortran -fPIC -w -O3 -shared -x f95 -o $(piclib * "." * Libdl.dlext) -`, "w") do f
        print(f, read(f90file, String))
    end
end

export f90_interpolation!

function f90_interpolation!( p :: Array{Float64,2}, fdtd :: FDTD )

    nx = Int32(fdtd.m.nx+1)
    ny = Int32(fdtd.m.ny+1)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(p)[2])

    ccall((:interpolation, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, fdtd.ebj, p)

end

export f90_deposition!

function f90_deposition!( fdtd :: FDTD, p :: Array{Float64, 2} )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = Int32(size(p)[2])

    f = fdtd.ebj

    ccall((:deposition, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx+1, ny+1, dx, dy, p, f)

    for i=1:nx+1
      fdtd.ebj[4:5,i,1]  .+= fdtd.ebj[4:5,i,ny+1]
      fdtd.ebj[4:5,i,ny+1]  .= fdtd.ebj[4:5,i,1]
    end
    for j=1:ny+1
      fdtd.ebj[4:5,1,j]  .+= fdtd.ebj[4:5,nx+1,j]
      fdtd.ebj[4:5,nx+1,j]  .= fdtd.ebj[4:5,1,j]
    end

    for i=1:nx, j=1:ny+1
       fdtd.jx[i,j] = 0.5 * (fdtd.ebj[4,i,j]+fdtd.ebj[4,i+1,j])
    end
    
    for i=1:nx+1, j=1:ny
       fdtd.jy[i,j] = 0.5 * (fdtd.ebj[5,i,j]+fdtd.ebj[5,i,j+1])
    end

end

export f90_push_x!

function f90_push_x!( p :: Array{Float64,2}, mesh :: Mesh, dt :: Float64)

    nbpart = Int32(size(p)[2])
    dimx = mesh.dimx
    dimy = mesh.dimy

    ccall((:push_x, piclib), Cvoid, (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}), nbpart, dimx, dimy, p, dt)

end

export f90_push_v!

function f90_push_v!( p :: Array{Float64, 2}, dt :: Float64)

    nbpart = Int32(size(p)[2])

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64}), nbpart, p, dt)

end
