module Equilibrium

using Interpolations
using Requires

const mu0 = 4*pi*1e-7

include("equil.jl")
export AxisymmetricEquilibrium
export Bfield, Efield, Jfield

include("limiter.jl")
export Limiter, in_vessel

@require EFIT begin
    include("efit_io.jl")
    export load, load_geqdsk, load_limiter, load_geqdsk_limiter
end

end # module
