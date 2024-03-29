using LinearAlgebra
using GEMPIC
using ParticleInCell
using Test

include("test_poisson.jl")
include("test_poisson_1d.jl")
include("test_landau_damping.jl")
include("test_push.jl")
include("test_fdtd.jl")
include("test_deposition.jl")
include("test_particles.jl")
include("test_particle_group.jl")
include("test_splitting.jl")
include("test_poisson_2d.jl")
include("test_particle_mesh_coupling_spline_2d.jl")
include("test_spline_pp.jl")
include("test_particle_mesh_coupling_spline_1d.jl")
include("test_maxwell_1d_fem.jl")
include("test_efd.jl")
include("test_particles.jl")
include("bupdate.jl")
