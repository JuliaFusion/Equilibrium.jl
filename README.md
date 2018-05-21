# Equilibrium.jl

Provides tools for working with MHD Equilibrium.

```julia
using EFIT
using Equilibrium

g = readg("g000001.01000")
M = AxisymmetricEquilibrium(g)
limiter = Limiter(g)

# or
# M, limiter = read_geqdsk("g000001.01000")

Bfield(M, r, z)

Jfield(M, r, z)

in_vessel(limiter,r,z)
```

