__precompile__()

module Equilibrium

using EFIT
using LinearAlgebra
using Interpolations
using ForwardDiff
using StaticArrays
using Trapz

import PolygonOps
import Contour
import Optim

const mu0 = 4*pi*1e-7

abstract type AbstractEquilibrium end
export AbstractEquilibrium

# Fallbacks
_not_implemented(M) = error("$(typeof(M)) has not implemented this functionality")
(M::AbstractEquilibrium)(x,y) = _not_implemented(M)
magnetic_axis(M::AbstractEquilibrium,r,z) = _not_implemented(M)
limits(M::AbstractEquilibrium) = _not_implemented(M)
psi_limits(M::AbstractEquilibrium) = _not_implemented(M)
psi_gradient(M::AbstractEquilibrium) = _not_implemented(M)
electric_potential(M::AbstractEquilibrium, psi) = _not_implemented(M)
electric_potential_gradient(M::AbstractEquilibrium) = _not_implemented(M)
pressure(M::AbstractEquilibrium, psi) = _not_implemented(M)
poloidal_current(M::AbstractEquilibrium, psi) = _not_implemented(M)
pressure_gradient(M::AbstractEquilibrium, psi) = _not_implemented(M)
poloidal_current_gradient(M::AbstractEquilibrium, psi) = _not_implemented(M)
cocos(M::AbstractEquilibrium) = _not_implemented(M)
B0Ip_sign(M::AbstractEquilibrium) = _not_implemented(M)

# Equilibrium API
export magnetic_axis, limits, psi_limits, psi_gradient, electric_potential, electric_potential_gradient
export pressure, poloidal_current, pressure_gradient, poloidal_current_gradient, safety_factor, cocos, B0Ip_sign
export plasma_current, beta_n

include("cocos.jl")
export COCOS, cocos, check_cocos, identify_cocos, transform_cocos
export cylindrical_cocos, poloidal_cocos, cylindrical_cocos_indices, poloidal_cocos_indices

include("boundary.jl")
export Boundary, PlasmaBoundary, FluxSurface, Wall, in_boundary, in_plasma, in_vessel
export boundary, flux_surface, circumference, average, area, area_average, volume, volume_average

include("shape.jl")
export PlasmaShape, MillerShape, MillerExtendedHarmonicShape, shape
export curvature, triangularity, squareness, tilt, elevation, ovality
export elongation, aspect_ratio, major_radius, minor_radius

include("solovev.jl")
export SolovevEquilibrium, solovev, clear_cache

include("fields.jl")
export Bfield, Efield, Jfield, EMFields, fields, gradB, curlB, poloidal_Bfield, poloidal_Jfield

include("efit.jl")
export EFITEquilibrium, efit

include("efit_io.jl")
export read_geqdsk

include("transp_io.jl")
export transp_potential!

end # module
