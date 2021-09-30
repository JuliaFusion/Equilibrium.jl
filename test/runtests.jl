using Test

using ForwardDiff
using Equilibrium
using EFIT

mu0 = 4pi*10^-7

@testset "Equilibrium" begin

include("solovev.jl")

include("efit.jl")

end
