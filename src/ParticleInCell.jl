module ParticleInCell

using OffsetArrays
import Sobol

include("simulation.jl")
include("particles.jl")
include("mesh.jl")
include("fields.jl")
include("decalage.jl")
include("maxwell.jl")
include("interpolate.jl")

end
