export Mesh

struct Mesh

    xmin :: Float64
    xmax :: Float64
    nx   :: Int32
    dx   :: Float64
    ymin :: Float64
    ymax :: Float64
    ny   :: Int32
    dy   :: Float64

    function Mesh( xmin, xmax, nx, ymin, ymax, ny )

        dx = (xmax - xmin) / nx
        dy = (ymax - ymin) / ny

        new( xmin, xmax, nx, dx, ymin, ymax, ny, dy )

    end

end

export MeshFields

mutable struct MeshFields

    mesh :: Mesh

    e :: Array{Float64,3}
    ρ :: Array{Float64,2}

    function MeshFields( mesh :: Mesh )
	    
	     nx, ny = mesh.nx, mesh.ny

         e = zeros(Float64, (2, nx+1, ny+1))
         ρ = zeros(Float64, (nx+1, ny+1))

         new( mesh, e, ρ)

    end

end

