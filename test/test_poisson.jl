@testset "Poisson 2D on rectangular grid" begin

    alpha = 0.05
    kx    = 0.5
    ky    = 1.0

    mesh = Mesh( 0, 2π/kx, 64, 0, 2π/ky, 128)

    fields = MeshFields( mesh )
    solutions = MeshFields( mesh)

    x = range(mesh.xmin, stop=mesh.xmax, length=mesh.nx+1) |> collect
    y = range(mesh.ymin, stop=mesh.ymax, length=mesh.ny+1) |> collect
    
    fields.ρ  .= - 8 * sin.(2*x) .* cos.(2*y')

    solutions.e[1,:,:] .=   2 * (cos.(2 .* x) .* cos.(2 .* y'))
    solutions.e[2,:,:] .= - 2 * (sin.(2 .* x) .* sin.(2 .* y'))

    poisson! = Poisson( mesh )

    poisson!( fields )

    err = errors( fields, solutions )

    println(" error : $err ")

    @test err ≈ 0.0 atol = 1e-14

    fields.ρ .= - 4 * ( sin.(2*x) .+ cos.(2*y') )

    poisson!( fields )

    nx, ny = mesh.nx, mesh.ny

    for j in 1:ny+1, i in 1:nx+1
        solutions.e[1,i,j] =   2*cos(2*x[i])
        solutions.e[2,i,j] = - 2*sin(2*y[j])
    end

    gnuplot("test2.dat", fields)
    gnuplot("solu2.dat", solutions)

    println(" error : $err ")
    err = errors( fields, solutions )

    @test err ≈ 0.0 atol = 1e-14

end
