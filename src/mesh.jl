export OneDGrid, TwoDGrid

"""
    TwoDGrid( dimx, nx, dimy, ny)

Generate a cartesians mesh on rectangle `dimx`x `dimy` with `nx` x `ny` points

- `nx` : indices are in [1:nx]
- `ny` : indices are in [1:ny]
- `dimx = xmax - xmin`
- `dimy = ymax - ymin`
- `x, y` : node positions
- `dx, dy` : step size
"""
struct TwoDGrid

    nx::Int
    ny::Int
    xmin::Float64
    xmax::Float64
    ymin::Float64
    ymax::Float64
    dimx::Float64
    dimy::Float64
    x::Vector{Float64}
    y::Vector{Float64}
    dx::Float64
    dy::Float64
    ex::Array{Float64,2}
    ey::Array{Float64,2}
    bz::Array{Float64,2}
    jx::Array{Float64,2}
    jy::Array{Float64,2}

    function TwoDGrid(xmin, xmax, nx, ymin, ymax, ny)

        dimx = xmax - xmin
        dimy = ymax - ymin

        x = LinRange(0, dimx, nx + 1) |> collect
        y = LinRange(0, dimy, ny + 1) |> collect

        dx = dimx / nx
        dy = dimy / ny

        ex = zeros(nx+1,ny+1)
        ey = zeros(nx+1,ny+1)
        bz = zeros(nx+1,ny+1)
        jx = zeros(nx+1,ny+1)
        jy = zeros(nx+1,ny+1)

        new(nx, ny, xmin, xmax, ymin, ymax, dimx, dimy, x, y, dx, dy, ex, ey, bz, jx, jy)

    end

end

TwoDGrid(dimx, nx :: Int, dimy, ny :: Int) = TwoDGrid(0.0, dimx, nx, 0.0, dimy, ny)

export OneDGrid

"""
    TwoDGrid( xmin, xmax, nx )

Simple structure to store mesh data from 1 to 3 dimensions
"""
struct OneDGrid 

    nx    :: Int
    xmin  :: Float64
    xmax  :: Float64
    Lx    :: Float64 
    dx    :: Float64

    function OneDGrid( xmin :: Real, xmax :: Real, nx :: Int) 

        Lx = xmax - xmin
        dx = Lx / (nx - 1)

        new( nx, xmin, xmax, Lx, dx )

    end

end 

export get_x

"""  
    get_x( mesh, i )

Get position
"""
function get_x( m :: OneDGrid, i )

    m.xmin + (i - 1) * m.dx
    
end 

"""  
    get_cell_and_offset( mesh, x )

Get cell and offset

We compute the cell indices where the particle is and its relative 
normalized position inside the cell

"""
function get_cell_and_offset( m :: OneDGrid, x ) :: Int64

    cell   = floor(Int64, ((x - m.xmin) / m.Lx) * m.nx) + 1
    offset = (x - get_x( m, cell)) / dx

	cell, offset
    
end 
