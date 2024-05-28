export MeshFields2D

struct MeshFields2D

    mesh :: TwoDGrid

    e :: Array{Float64,3}
    ρ :: Array{Float64,2}

    function MeshFields2D( mesh :: TwoDGrid )
	    
	     nx, ny = mesh.nx, mesh.ny

         e = zeros(Float64, (2, nx+1, ny+1))
         ρ = zeros(Float64, (nx+1, ny+1))

         new( mesh, e, ρ)

    end

end

using FFTW

export Poisson2D

struct Poisson2D

    mesh   :: TwoDGrid
    kx     :: Array{Float64, 2}
    ky     :: Array{Float64, 2}
    ρ̃      :: Array{ComplexF64, 2}
    fft_ex :: Array{ComplexF64, 2}
    fft_ey :: Array{ComplexF64, 2}

    function Poisson2D( mesh :: TwoDGrid )

        nx = mesh.nx
        ny = mesh.ny

        kx0 = 2π / (mesh.xmax - mesh.xmin)
        ky0 = 2π / (mesh.ymax - mesh.ymin)

        kx = zeros(Float64, (nx÷2+1,ny))
        ky = zeros(Float64, (nx÷2+1,ny))

        for ik=1:nx÷2+1
           kx1 = (ik-1)*kx0
           for jk = 1:ny÷2
              kx[ik,jk] = kx1
              ky[ik,jk] = (jk-1)*ky0
           end
           for jk = ny÷2+1:ny
              kx[ik,jk] = kx1
              ky[ik,jk] = (jk-1-ny)*ky0
           end 
        end

        kx[1,1] = 1.0
        k2  = kx .* kx .+ ky .* ky
        kx .= kx ./ k2
        ky .= ky ./ k2

        ρ̃ = zeros(ComplexF64,(nx÷2+1,ny))

        fft_ex = similar(ρ̃)
        fft_ey = similar(ρ̃)

        new( mesh, kx, ky, ρ̃, fft_ex, fft_ey)

    end	    

end


"""

   poisson!( fields )

Solve the equation Δ Φ = - fields.ρ

 fields.ex = ∂ Φ / ∂ x
 fields.ey = ∂ Φ / ∂ y 


"""
function solve!( fields :: MeshFields2D, p :: Poisson2D )

    nx, ny = p.mesh.nx, p.mesh.ny
    dx, dy = p.mesh.dx, p.mesh.dy

    p.ρ̃ .= rfft(view(fields.ρ,1:nx,1:ny))

    p.fft_ex .= -1im .* p.kx .* p.ρ̃
    p.fft_ey .= -1im .* p.ky .* p.ρ̃

    fields.e[1,1:nx,1:ny] .= irfft(p.fft_ex, nx)
    fields.e[2,1:nx,1:ny] .= irfft(p.fft_ey, nx)

    fields.e[1,nx+1,:] .= view(fields.e,1,1,:)
    fields.e[1,:,ny+1] .= view(fields.e,1,:,1)
    fields.e[2,nx+1,:] .= view(fields.e,2,1,:)
    fields.e[2,:,ny+1] .= view(fields.e,2,:,1)

    @views sum(  fields.e[1,:,:] .* fields.e[1,:,:] 
              .+ fields.e[2,:,:] .* fields.e[2,:,:]) * dx * dy

end
