# ParticleInCell.jl

[![Documentation](https://github.com/juliavlasov/ParticleInCell.jl/workflows/Documentation/badge.svg)](https://juliavlasov.github.io/ParticleInCell.jl/dev)

Particle In Cell code in Julia

Work in progress, you can find first tests in the documentation.

I am working on this Julia code to compare performances with a Fortran code that 
solves same equation with same parameters and same numerical method.
You can find it [here](https://github.com/pnavaro/vm_nonunif). This is an old code written in 2005
and not well optimized but it takes 6 seconds with 204800 particles.
I will note here times of this Julia code and what I have done to speed-up things.

- 323 seconds : First version 
- 303 seconds : Use julia -O3 --check-bounds=no
- 302 seconds : change shape of positions and velocities arrays for particles
- 156 seconds : Regroup ex, ey, and bz in a same array eb(3,nx,ny)

## Other Julia PIC codes 

- [ParticleInCell.jl](https://github.com/adamslc/ParticleInCell.jl) by Luke Adams.
- [Julia discourse thread](https://discourse.julialang.org/t/pic-particle-in-cell-space-charge-tracking-simulation/)
