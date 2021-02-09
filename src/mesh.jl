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

    nx::Any
    ny::Any
    dimx::Any
    dimy::Any
    x::Any
    y::Any
    dx::Any
    dy::Any

    function Mesh(dimx, nx, dimy, ny)

        x = LinRange(0, dimx, nx + 1) |> collect
        y = LinRange(0, dimy, ny + 1) |> collect

        dx = dimx / nx
        dy = dimy / ny

        new(nx, ny, dimx, dimy, x, y, dx, dy)

    end

end
