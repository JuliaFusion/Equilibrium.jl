__precompile__()

module Equilibrium

using EFIT
using LinearAlgebra
using Interpolations
using StaticArrays

const mu0 = 4*pi*1e-7

include("equil.jl")
export AxisymmetricEquilibrium
export Bfield, Efield, Jfield, EMFields, fields

include("limiter.jl")
export Limiter, in_vessel

include("efit_io.jl")
export read_geqdsk

include("transp_io.jl")
export transp_potential!

end # module
