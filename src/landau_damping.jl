using Sobol

export LandauDamping

"""
    Landau( α, kx)

Test structure to initialize a particles distribtion for
Landau damping test case in 1D1V and 1D2V

"""
struct LandauDamping

    alpha::Float64
    kx::Float64

end

function newton(r, alpha, kx)
    x0, x1 = 0.0, 1.0
    r *= 2π / kx
    while (abs(x1 - x0) > 1e-12)
        p = x0 + alpha * sin(kx * x0) / kx
        f = 1 + alpha * cos(kx * x0)
        x0, x1 = x1, x0 - (p - r) / f
    end
    x1
end

export sample!


"""
    sample!(d, pg)

Sampling from a probability distribution to initialize a Landau damping in
1D2V space.

```math
f_0(x,v,t) = \\frac{n_0}{2π v_{th}^2} ( 1 + \\alpha cos(k_x x))
 exp( - \\frac{v_x^2+v_y^2}{2 v_{th}^2})
```
The newton function solves the equation ``P(x)-r=0`` with Newton’s method
```math
x^{n+1} = x^n – (P(x)-(2\\pi r / k)/f(x) 
```
with 
```math
P(x) = \\int_0^x (1 + \\alpha cos(k_x y)) dy = x + \\frac{\\alpha}{k_x} sin(k_x x)
```
"""
function sample!(pg::ParticleGroup{1,2}, d::LandauDamping)

    alpha, kx = d.alpha, d.kx

    s = SobolSeq(1)
    nbpart = pg.n_particles

    for i = 1:nbpart
        v = sqrt(-2 * log((i - 0.5) / nbpart))
        x = Sobol.next!(s)
        set_x!(pg, i, newton(x, alpha, kx))
        set_v!(pg, i, v)
        set_w!(pg, i, 2 * pi / kx / nbpart)
    end

end

"""
    sample!(d, pg)

Sampling from a probability distribution to initialize a Landau damping in
2D2V space.

```math
f_0(x,v,t) = \\frac{n_0}{2π v_{th}^2} ( 1 + \\alpha cos(k_x x))
 exp( - \\frac{v_x^2+v_y^2}{2 v_{th}^2})
```
"""
function sample!(pg::ParticleGroup{2,2}, mesh :: TwoDGrid, d::LandauDamping)


    alpha, kx = d.alpha, d.kx

    @assert mesh.dimx ≈ 2π / kx

    s = Sobol.SobolSeq(3)
    nbpart = pg.n_particles

    for i = 1:nbpart
        v = sqrt(-2 * log((i - 0.5) / nbpart))
        x, y, θ = Sobol.next!(s)
        θ = θ * 2π
        pg.array[1, i] = newton(x, alpha, kx)
        pg.array[2, i] = mesh.ymin + y * mesh.dimy
        pg.array[3, i] = v * cos(θ)
        pg.array[4, i] = v * sin(θ)
        pg.array[5, i] = (mesh.dimx * mesh.dimy) / nbpart
    end

end
