export MeshFields

struct MeshFields

    mesh :: TwoDGrid

    e :: Array{Float64,3}
    ρ :: Array{Float64,2}

    function MeshFields( mesh :: TwoDGrid )
	    
	     nx, ny = mesh.nx, mesh.ny

         e = zeros(Float64, (2, nx+1, ny+1))
         ρ = zeros(Float64, (nx+1, ny+1))

         new( mesh, e, ρ)

    end

end

