# ParticleInCell.jl

[![Documentation](https://github.com/juliavlasov/ParticleInCell.jl/workflows/Documentation/badge.svg)](https://juliavlasov.github.io/ParticleInCell.jl/dev)

Particle In Cell code in Julia

Work in progress, you can find first tests in the documentation.

I am working on this Julia code to compare performances with a Fortran code that 
solves same equation with same parameters and same numerical method.
You can find it [here](https://github.com/pnavaro/vm_nonunif). This is an old code written in 2005
and not well optimized but it takes 6 seconds with 204800 particles.
I note here times of this Julia code and what I have done to speed-up things.

- 323 seconds : First version 
- 303 seconds : Use julia -O3 --check-bounds=no
- 302 seconds : change shape of positions and velocities arrays for particles
- 156 seconds : Regroup ex, ey, and bz in a same array eb(3,nx,ny)
- 176 seconds : put the array ebj in fdtd type
- 140 seconds : put the particles data in one array, the particle push_v is now very fast.
- 077 seconds : replace the julia interpolation function by a call to the fortran subroutine.
- 019 seconds : replace the julia deposition by a call to the fortran subroutine
- 006 seconds : vectorize and use views in function push_x instead of a loop
- 135 seconds : back to julia functions for particles but without `mod1` calls.
- 003 seconds : back to julia functions, I set the types of struct members.

Now I increase the number of particles to 1024000 and begin to parallelize particles motion.

- 58 seconds : First serial time, It spends 24 seconds in `push_x`, we can do better.
- 46 seconds : Move interpolation step inside the `push_v` function and reduce memory print. 
- 31 seconds : Add @threads in `push_x` and `push_v` functions. 
- 17 seconds : by adding @threads in `compute_current`

## Other Julia PIC codes 

- [ParticleInCell.jl](https://github.com/adamslc/ParticleInCell.jl) by Luke Adams.
- [Julia discourse thread](https://discourse.julialang.org/t/pic-particle-in-cell-space-charge-tracking-simulation/)
