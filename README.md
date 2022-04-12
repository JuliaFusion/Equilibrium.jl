# Equilibrium.jl

Equilibrium.jl provides tools for working with solutions of the Grad-Shafranov Equation.

## AbstractEquilibrium API

```julia
using Equilibrium
typeof(S) <: AbstractEquilibrium

psi = S(r, z) # Poloidal flux at r,z
gradpsi = psi_gradient(S, r, z)

B = Bfield(S, r, z)
Bp = poloidal_Bfield(S, r, z)

J = Jfield(S, r, z)
Jp = poloidal_Jfield(S, r, z)

F = poloidal_current(S, psi)
Fprime = poloidal_current_gradient(S, psi)

p = pressure(S, psi)
pprime = pressure_gradient(S, psi)

V = electric_potential(S, psi)
gradV = electric_potential_gradient(S, psi)

q = safety_factor(S, psi)

maxis = magnetic_axis(S)

btip = B0Ip_sign(S)

rlims, zlims = limits(S)

psi_lims = psi_limits(S)

cc = cocos(S) # Return COCOS structure

fs = flux_surface(S, psi) # returns a boundary object

```

## Solov'ev Equilibrium
Solov'ev Equilibrium are analytic solutions to the Grad-Shafranov equation where the p' and FF' are constant.
The resulting Grad-Shafranov equation takes the form `Δ⋆ψ = α + (1-α)x²` where `α` is some constant.
The boundary conditions are found using a plasma shape parameterization.

```julia

# ITER parameters
δ = 0.33        # Triangularity
ϵ = 0.32        # Inverse aspect ratio a/R0
κ = 1.7         # Elongation
B0 = 5.3        # Magnitude of Toroidal field at R0 [T]
R0 = 6.2        # Major Radius [m]
qstar = 1.57    # Kink safety factor
alpha = -0.155  # constant

S = solovev(B0, MillerShape(R0, ϵ, κ, δ), alpha, qstar, B0_dir=1, Ip_dir=1)


SolovevEquilibrium
  B0 = 2.0 [T]
  S  = MillerShape{Float64}(6.2, 0.32, 1.7, 0.33)
  α  = -0.155
  q⋆ = 1.57
  βp = 1.1837605469381924
  βt = 0.049177281028224634
  σ  = 1
  diverted  = false
  symmetric = true
```

## EFIT Equilibrium
EFIT geqdsk files are a commonly used file format.
Here we provide routines for converting the GEQDSK files into an Equilibrium object.

```julia
using EFIT

g = readg("g000001.01000")
M = efit(g, clockwise_phi=false) # direction of phi needed to determine COCOS ID
wall = Wall(g)

in_vessel(wall, r, z)

# or
# M, wall = read_geqdsk("g000001.01000",clockwise_phi=false)

```

## COCOS: Tokamak Coordinate Conventions
We provide routines for working determining, transforming, and checking COCOS's.
```julia
julia> cocos(3)
COCOS = 3
 e_Bp  = 0
 σ_Bp  = -1
 σ_RΦZ = (R,Φ,Z): 1
 σ_ρθΦ = (ρ,Φ,θ): -1
 Φ from top: CCW
 θ from front: CCW
 ψ_ref: Decreasing assuming +Ip, +B0
 sign(q) = -1 assuming +Ip, +B0
 sign(p') = 1 assuming +Ip, +B0

julia> transform_cocos(3,1)
Dict{Any, Any} with 14 entries:
  "Z"        => 1.0
  "Q"        => -1
  "P"        => 1.0
  "B"        => 1.0
  "F_FPRIME" => -1.0
  "ψ"        => -1.0
  "TOR"      => 1.0
  "Φ"        => 1.0
  "PSI"      => -1.0
  "I"        => 1.0
  "J"        => 1.0
  "R"        => 1.0
  "F"        => 1.0
  "PPRIME"   => -1.0
```

## Boundaries
Equilibrium.jl also provides routines for working with boundries such as walls or flux surfaces. Internally boundaries are stored as a list of points forming a polygon.

```julia

fs = flux_surface(S, psi) # alternatively fs = boundary(S, psi)

in_plasma(fs, r, z) # or in_vessel(fs, r, z), in_boundary(fs, r, z)

cicumference(fs)

area(fs) # Area enclosed by the boundary

volume(fs) # assuming toroidal symmetry. F can be a vector with the same length as fs or a function of (r,z)

average(fs, F) # Average F over the boundary

area_average(fs, F) # average F over the area

volume_average(fs, F) # average F over the volume
```
