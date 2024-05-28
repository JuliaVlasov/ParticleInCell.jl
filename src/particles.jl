export Particle, Particles

struct Particle

    x :: Float64
    v :: ComplexF64
    e :: ComplexF64
    b :: Float64

end

mutable struct Particles

    nbpart :: Int64

    x :: Array{Float64,2}
    v :: Array{Float64,2}
    e :: Array{Float64,2}
    b :: Vector{Float64}
    t :: Vector{Float64}
    w :: Float64

    function Particles( nbpart :: Int64, w :: Float64 )

        x = zeros(Float64, (2,nbpart))
        v = zeros(Float64, (2,nbpart))
        e = zeros(Float64, (2,nbpart))
        b = zeros(Float64, nbpart)
        t = zeros(Float64, nbpart)

        new( nbpart, x, v, e, b, t , w )

    end

end

export landau_sampling

function landau_sampling( mesh :: TwoDGrid, alpha, kx, nbpart :: Int64 )

    nx, ny = mesh.nx, mesh.ny
    dx, dy = mesh.dx, mesh.dy
    dimx = mesh.xmax - mesh.xmin
    dimy = mesh.ymax - mesh.ymin

    weight = (dimx * dimy) / nbpart

    particles = Particles(nbpart, weight )

    function newton(r)
        x0, x1 = 0.0, 1.0
        r *= 2π / kx
        while (abs(x1-x0) > 1e-12)
            p = x0 + alpha * sin( kx * x0) / kx 
            f = 1 + alpha * cos( kx * x0)
            x0, x1 = x1, x0 - (p - r) / f
        end
        x1
    end
    
    s = Sobol.SobolSeq(3)

    for i=1:nbpart

        v = sqrt(-2 * log( (i-0.5)/nbpart))
        r1, r2, r3  = Sobol.next!(s)
        θ = r1 * 2π
        particles.x[1,i] = mesh.xmin + newton(r2) * dimx
        particles.x[2,i] = mesh.ymin + r3 * dimy
        particles.v[1,i] = v * cos(θ)
        particles.v[2,i] = v * sin(θ)

    end

    particles

end


