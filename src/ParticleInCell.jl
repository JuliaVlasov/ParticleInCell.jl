module ParticleInCell

struct Particle

   cell :: Int32
   dx   :: Float32
   v    :: Float64
   w    :: Float32

end

include("simulation.jl")
