# Equilibrium.jl

Provides tools for working with MHD Equilibrium.

```julia
using EFIT
using Equilibrium

g = readg("g000001.01000")

M = load(g) # or M = load_geqdsk("g000001.01000")

limiter = load_limiter(g) # or load_geqdsk_limiter("g000001.01000")

Bfield(M, r, z)

Jfield(M, r, z)

in_vessel(limiter,r,z)
```

