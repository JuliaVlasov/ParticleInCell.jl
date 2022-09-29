export gnuplot
export errors

function gnuplot( filename :: String, fields :: MeshFields )

    open(filename, "w") do f

        nx = fields.mesh.nx
        ny = fields.mesh.ny
        dx = fields.mesh.dx
        dy = fields.mesh.dy
    
        for i in 1:nx+1
            for j in 1:ny+1
                write(f, 
	    	        string((i-1)*dx),        "  ", 
	    	        string((j-1)*dy),        "  ",  
	    	        string(fields.e[1,i,j]), "  ", 
	    	        string(fields.e[2,i,j]), "  ",
	    	        string(fields.ρ[i,j]),  "\n")
            end
	        write(f,"\n")
        end

    end

end

function errors( computed :: MeshFields, reference :: MeshFields )

    err = maximum(abs.(computed.e .- reference.e))

    err

end

