module ParticleInCell

using LinearAlgebra
using Random
using SparseArrays
import Sobol
import GEMPIC: OneDGrid, TwoDGrid, ParticleGroup, ParticleMeshCoupling1D, ParticleMeshCoupling2D
export OneDGrid, TwoDGrid, ParticleGroup, ParticleMeshCoupling1D, ParticleMeshCoupling2D
import GEMPIC: add_charge!
export add_charge!

include("particles.jl")
include("particle_1d1v.jl")
include("landau_damping.jl")
include("poisson_1d.jl")
include("poisson_2d.jl")
include("cloud_in_cell.jl")
include("pushers.jl")
include("pushers_m6.jl")
include("fdtd.jl")
include("pic.jl")
include("ua_type.jl")
include("meshfields.jl")
include("read_particles.jl")
include("plasma.jl")
include("integrate.jl")
include("gnuplot.jl")
include("ua_steps.jl")
include("compute_rho_cic.jl")
include("interpolate_cic.jl")
include("compute_rho_m6.jl")
include("interpolation_m6.jl")

end
