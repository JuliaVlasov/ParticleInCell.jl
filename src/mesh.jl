export Mesh

"""
    Mesh( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [0:nx]
- `ny` : indices are in [0:ny]
- `dimx` x `dimy`: mesh area
- `x, y` : node positions
- `dx, dy` : step size
"""
struct Mesh

    nx
    ny
    dimx
    dimy
    x
    y
    dx
    dy

    function Mesh(dimx, nx, dimy, ny)

        x = LinRange(0, dimx, nx+1) |> collect
        y = LinRange(0, dimy, ny+1) |> collect
        
        dx = dimx / nx
        dy = dimy / ny
        
        new( nx, ny, dimx, dimy, x, y, dx, dy )

    end

end
