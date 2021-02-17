export Mesh

"""
    Mesh( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [1:nx]
- `ny` : indices are in [1:ny]
- `dimx` x `dimy`: mesh area
- `x, y` : node positions
- `dx, dy` : step size
"""
struct Mesh

    nx::Int
    ny::Int
    dimx::Float64
    dimy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    dx::Float64
    dy::Float64
    ex::Array{Float64,2}
    ey::Array{Float64,2}
    bz::Array{Float64,2}

    function Mesh(dimx, nx, dimy, ny)

        x = LinRange(0, dimx, nx + 1) |> collect
        y = LinRange(0, dimy, ny + 1) |> collect

        dx = dimx / nx
        dy = dimy / ny

        ex = zeros(nx+1,ny+1)
        ey = zeros(nx+1,ny+1)
        bz = zeros(nx+1,ny+1)

        new(nx, ny, dimx, dimy, x, y, dx, dy, ex, ey, bz)

    end

end
