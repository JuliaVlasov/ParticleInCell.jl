using FFTW

export Maxwell

struct Maxwell

    kx::Vector{Float64}
    ky::Vector{Float64}

    function Maxwell(mesh)

        nx = mesh.nx
        kx = 2π / mesh.dimx .* vcat(0:nx÷2-1, -nx÷2:-1)
        ny = mesh.ny
        ky = 2π / mesh.dimy .* vcat(0:ny÷2-1, -ny÷2:-1)

        new(kx, ky)

    end

end

export faraday!, ampere_maxwell!

"""
    ampere_maxwell!( ex, ey, solver, bz, dt)

```math
Ex^{t+dt} = Ex^{t} + dt \\big( \\frac{\\partial Bz}{\\partial y} - Jx \\big)
```
    
```math
Ey^{t+dt} = Ey^{t} - dt \\big( \\frac{\\partial Bz}{\\partial x} - Jy \\big)
```
"""
function ampere_maxwell!(ex, ey, s::Maxwell, bz, jx, jy, dt)

    dbz_dy = real(ifft(1im .* s.ky' .* fft(bz, 2)))
    ex .= ex .+ dt .* dbz_dy .- dt .* jx

    dbz_dx = real(ifft(1im .* s.kx .* fft(bz, 1)))
    ey .= ey .- dt .* dbz_dx .- dt .* jy

end


"""
    faraday!( bz, solver, ex, ey, dt)

```math
Bz^{t+dt} = Bz^{t} + dt (  \\frac{\\partial Ey}{\\partial x}
-   \\frac{\\partial Ex}{\\partial y})
"""
function faraday!(bz, s::Maxwell, ex, ey, dt::Float64)

    dex_dy = real(ifft(1im .* s.ky' .* fft(ex, 2)))
    dey_dx = real(ifft(1im .* s.kx .* fft(ey, 1)))

    bz .+= dt .* (dex_dy .- dey_dx)

end
