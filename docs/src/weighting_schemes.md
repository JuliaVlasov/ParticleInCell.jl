# Weighting schemes

This a comparison between four weighting schemes using the Landau damping test case.

- NGP : Nearest Grid Point
- CIC : Clud In Cell
- TSC : Triangular Shape Cloud
- M4 : Quartic Spline kernel (Monaghan)
- M6 : Quintic spline kernel (Monaghan)

```@example ws
using ParticleInCell
using Plots
using DispersionRelations
```

```@example ws
nx = 128
ny = 8
alpha = 0.1
kx = 0.5
ky = 0.0
dimx = 2pi / kx
dimy = 1.0
dt = 0.1
```

```@example ws
function simulation( kernel)

    nbpart = 100 * nx * ny
    mesh = TwoDGrid(dimx, nx, dimy, ny)

    p = ParticleGroup{2,2}(nbpart, charge = 1.0, mass = 1.0, n_weights = 1)

    sampler = LandauDamping(alpha, kx)
    ParticleInCell.sample!(p, mesh, sampler)
    
    ρ_ref = zeros(nx, ny)

    for i = 1:nx, j = 1:ny
        ρ_ref[i, j] = alpha * cos(kx * mesh.x[i])
    end
    
    ρ = compute_rho(p, kernel, mesh)
    
    @show sum((ρ .- ρ_ref).^2)
    
    ex = similar(ρ)
    ey = similar(ρ)
    bz = zero(ρ)
    
    poisson = TwoDPoissonPeriodic(mesh)
    
    solve!(ex, ey, poisson, ρ)
    
    nsteps = 200
    
    energy(ex) = sqrt(sum( ex .^ 2) * mesh.dx * mesh.dy)
    
    t = [0.0]
    e = [energy(ex)]
    
    for i in 1:nsteps
        
        push_x!(p, mesh, 0.5dt)
        
        ρ .= compute_rho(p, kernel, mesh)
        solve!(ex, ey, poisson, ρ)
        push!(t, i*dt)
        push!(e, energy(ex))
        
        push_v!(p, kernel, mesh, ex, ey, bz, dt)
    
        push_x!(p, mesh, 0.5dt)
        
    end
    
    t, e
    
end
```

```@example ws
@time t, e = simulation(NearestGridPoint())
line, γ = fit_complex_frequency(t, e)
plot(t, e, yscale = :log10, label="NGP")  
plot!(t, line, label = "$(imag(γ))")
@time t, e = simulation(CloudInCell())
line, γ = fit_complex_frequency(t, e)
plot!(t, e, yscale = :log10, label="CIC")  
plot!(t, line, label = "$(imag(γ))")
@time t, e = simulation(TriangularShapeCloud())
line, γ = fit_complex_frequency(t, e)
plot!(t, e, yscale = :log10, label="TSC")  
plot!(t, line, label = "$(imag(γ))")
@time t, e = simulation(M4())
line, γ = fit_complex_frequency(t, e)
plot!(t, e, yscale = :log10, label="M4")  
plot!(t, line, label = "$(imag(γ))")
@time t, e = simulation(M6())
line, γ = fit_complex_frequency(t, e)
plot!(t, e, yscale = :log10, label="M6")  
plot!(t, line, label = "$(imag(γ))")
```
