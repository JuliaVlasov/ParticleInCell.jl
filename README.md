# ParticleInCell.jl

[![Documentation](https://github.com/juliavlasov/ParticleInCell.jl/workflows/Documentation/badge.svg)](https://juliavlasov.github.io/ParticleInCell.jl/dev)

Particle In Cell code in Julia

This software is developed in order to experiment optimizations and to try the best way to have an efficient package for particle-mesh numerical methods for plasmas.

It is a work in progress, you can find first tests in the documentation.

I compared performances of this Julia code with a Fortran version that 
solves same equation with same parameters and same numerical method.
You can find it [here](https://github.com/pnavaro/vm_nonunif). This is an old code written in 2005 with Régine Barthelmé and even is is not well optimized, it takes only 6 seconds with 204800 particles. The Julia version takes 3 seconds.

## Other Julia PIC codes 

- [ParticleInCell.jl](https://github.com/JuliaPlasma/ParticleInCell.jl) by Luke Adams.
- [MPIParticleInCell.jl](https://github.com/adamslc/MPIParticleInCell.jl) by Luke Adams.
- [ElectrostaticPIC1D.jl](https://github.com/jwscook/ElectrostaticPIC1D.jl) by James Cook.
- [Julia discourse thread](https://discourse.julialang.org/t/pic-particle-in-cell-space-charge-tracking-simulation/)


# Papers

The simulations in the following papers have been made with this code. Some of them are implemented in FORTRAN translated in Julia. Both codes are in the repository. The interface and code are messy, sorry for that. Don't hésitate to contact me or use the [Discussions tab](https://github.com/JuliaVlasov/ParticleInCell.jl/discussions).

- P. Chartier, N. Crouseilles, M. Lemou, F. Méhats, X. Zhao, Uniformly accurate methods for three dimensional Vlasov equations under strong magnetic field with varying direction, SIAM J. Sc. Comput. 42, pp. 520-547, 2020.
- P. Chartier, N. Crouseilles, M. Lemou, F. Méhats, X. Zhao, Uniformly accurate methods for Vlasov equations with non-homogeneous strong magnetic field, Math. Comp. 88, pp. 2697-2736, 2019.
- P. Chartier, N. Crouseilles, X. Zhao, Numerical methods for the two-dimensional Vlasov-Poisson models in the finite Larmor radius approximation regime, J. Comput. Phys., 375, pp. 619-640, 2018.
- N. Crouseilles, S. Hirstoaga, X. Zhao, Multiscale Particle-in-Cell methods and comparisons for the long time two-dimensional Vlasov-Poisson equation with strong magnetic field, Comput. Phys. Comm. 222, pp. 136-151, 2018.
- N. Crouseilles, M. Lemou, F. Méhats, X. Zhao, Uniformly accurate Particle-In-Cell method for the long time solution of the two-dimensional Vlasov-Poisson equation with uniform strong magnetic field, J. Comput. Phys. 346, pp. 172-190, 2017.
- N. Crouseilles, M. Lemou, F. Méhats, X. Zhao, Uniformly accurate forward semi-Lagrangian methods for highly oscillatory Vlasov-Poisson equations, SIAM MMS 15, pp. 723-744, 2017.
