using Equilibrium
using Test
using EFIT

mu0 = 4pi*10^-7

@testset "Equilibrium" begin

include("numerical.jl")

include("solovev.jl")

end
