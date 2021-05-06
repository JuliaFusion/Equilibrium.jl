__precompile__()

module Equilibrium

using EFIT
using LinearAlgebra
using Interpolations
using ForwardDiff
using StaticArrays
using PolygonOps

const mu0 = 4*pi*1e-7

include("cocos.jl")
export COCOS, cocos, check_cocos, identify_cocos, transform_cocos
export cylindrical_cocos, poloidal_cocos, cylindrical_cocos_indices, poloidal_cocos_indices

include("equil.jl")
export AxisymmetricEquilibrium
export Bfield, Efield, Jfield, EMFields, fields, gradB, curlB

include("wall.jl")
export Wall, Wall2D, in_vessel

include("boundary.jl")
export PlasmaBoundary, in_plasma, circumference, area, area_average, volume, volume_average

include("efit_io.jl")
export read_geqdsk

include("transp_io.jl")
export transp_potential!

end # module
