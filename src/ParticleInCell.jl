module ParticleInCell

using LinearAlgebra
using Random
using SparseArrays
using OffsetArrays
import Sobol

include("particle_1d1v.jl")
include("particles.jl")
include("mesh.jl")
include("compute_rho.jl")
include("fdtd.jl")
include("pic.jl")

end
