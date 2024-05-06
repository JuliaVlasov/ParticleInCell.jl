# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,jl:light
#     text_representation:
#       extension: .jl
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Julia 1.10.3
#     language: julia
#     name: julia-1.10
# ---

# + endofcell="--"
using Sobol
using Plots

function test( nbpart :: Int64)
    
   up = Float64[]
   vp = Float64[]
    
   s = SobolSeq(1)

   for k=0:nbpart-1

      v = sqrt(-2 * log( (k+0.5)/nbpart))
      r = next!(s)
      θ = r[1] * 2π
      push!(up,  v * cos(θ))
      push!(vp,  v * sin(θ))

   end

   up, vp

end

@time up, vp = test(1000000);
# -

p = plot(layout=(2,1))
histogram!(p[1], up, normalize=true, bins=100)
plot!(p[1], x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
histogram!(p[2], vp, normalize=true, bins=100)
plot!(p[2], x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
# --

# + endofcell="--"
using Sobol
using Plots

function test_with_weights( nbpart :: Int64)
    
   xp = Float64[]
   wp = Float64[]
    
   s = SobolSeq(1)

   for k=0:nbpart-1

      r = next!(s)
      x = -6 + r[1] * 12.0
      push!(xp,  x)
      push!(wp, exp( - 0.5 * (x^2) ) / sqrt(2π))

   end

   xp, wp

end

@time xp, wp = test_with_weights(1000000);
# -

histogram(xp, weights=wp, normalize=true, bins=100)
plot!(x-> (exp(-x^2/2))/sqrt(2π), -6, 6)
# --


