module ParticleInCell

using LinearAlgebra
using Random
using SparseArrays
import Sobol
import GEMPIC: OneDGrid, TwoDGrid, ParticleGroup, ParticleMeshCoupling1D, ParticleMeshCoupling2D
export OneDGrid, TwoDGrid, ParticleGroup, ParticleMeshCoupling1D, ParticleMeshCoupling2D

include("particles.jl")
include("particle_1d1v.jl")
include("landau_damping.jl")
include("poisson_1d.jl")
include("poisson_2d.jl")
include("particle_mesh.jl")
include("pushers.jl")
include("fdtd.jl")
include("pic.jl")
include("ua_type.jl")
include("meshfields.jl")
include("read_particles.jl")
include("plasma.jl")
include("compute_rho.jl")
include("interpolation.jl")
include("integrate.jl")
include("gnuplot.jl")
include("poisson.jl")
include("ua_steps.jl")
include("landau.jl")
include("compute_rho_cic.jl")

end
