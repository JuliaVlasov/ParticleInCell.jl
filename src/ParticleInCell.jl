module ParticleInCell

using LinearAlgebra
using Random
using SparseArrays
import Sobol

include("mesh.jl")
include("particle_1d1v.jl")
include("particles.jl")
include("compute_rho.jl")
include("fdtd.jl")
include("particle_mesh.jl")
include("pic.jl")
include("low_level_bsplines.jl")
include("maxwell_1d_fem.jl")
include("splinepp.jl")
include("particle_group.jl")
include("particle_mesh_coupling_1d.jl")
include("particle_mesh_coupling_2d.jl")

end
