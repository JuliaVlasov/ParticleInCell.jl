const piclib = joinpath(@__DIR__,"fortran", "libpic.dylib")

export vm2d2v

function vm2d2v( nstep :: Int, time :: Vector{Float64}, energy :: Vector{Float64} )

    nstep = Int32(nstep)

    ccall((:vm2d2v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ptr{Float64}), nstep, time, energy)

end

export interpolation!

function interpolation!( particles :: Particles, fdtd :: FDTD )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = particles.nbpart

    f = fdtd.ebj
    p = particles.data

    ccall((:interpolation, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, f, p)

end

export deposition!

function deposition!( fdtd :: FDTD, particles :: Particles )

    nx = Int32(fdtd.m.nx)
    ny = Int32(fdtd.m.ny)
    dx = fdtd.m.dx
    dy = fdtd.m.dy
    nbpart = particles.nbpart

    f = fdtd.ebj
    p = particles.data
    jx = fdtd.jx
    jy = fdtd.jy

    ccall((:deposition, piclib), Cvoid, (Ref{Int32}, Ref{Int32}, Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}), nbpart, nx, ny, dx, dy, p, f, jx, jy)

end

function push_x!( p :: Particles, dimx, dimy, dt)

    nbpart = Int32(p.nbpart)

    ccall((:push_x, piclib), Cvoid, (Ref{Int32}, Ref{Float64}, Ref{Float64}, Ptr{Float64}, Ref{Float64}), nbpart, dimx, dimy, p.data, dt)

end

function push_v!( p :: Particles, nbpart, dt :: Float64)

    nbpart = Int32(nbpart)

    ccall((:push_v, piclib), Cvoid, (Ref{Int32}, Ptr{Float64}, Ref{Float64}), nbpart, p.data, dt)

end
