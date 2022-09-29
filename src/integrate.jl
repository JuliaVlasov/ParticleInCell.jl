export integrate

function integrate( field :: Array{Float64,2}, mesh :: Mesh )

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy

    sum(view(field,1:nx,1:ny)) * dx * dy

end
