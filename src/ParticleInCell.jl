module ParticleInCell

using LinearAlgebra
using Random
using SparseArrays
import Sobol
import GEMPIC: OneDGrid, TwoDGrid, ParticleGroup
import GEMPIC: ParticleMeshCoupling1D
import GEMPIC: ParticleMeshCoupling2D

include("particle_1d1v.jl")
include("landau_damping.jl")
include("poisson_1d.jl")
include("poisson_2d.jl")
include("particle_mesh.jl")
include("particles.jl")
include("fdtd.jl")
include("pic.jl")

end
