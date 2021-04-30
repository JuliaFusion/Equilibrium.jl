using Equilibrium
using Test
using EFIT

dir = @__DIR__
g1 = readg(dir*"/g150219.03200") # (+Bt, +Ip)
g2 = readg(dir*"/g133221.01150") # (-Bt, +Ip)
g3 = readg(dir*"/g159177.02700") # (+Bt, -Ip)
g4 = readg(dir*"/g153298.04400") # (-Bt, -Ip)

M1 = AxisymmetricEquilibrium(g1)
M2 = AxisymmetricEquilibrium(g2)
M3 = AxisymmetricEquilibrium(g3)
M4 = AxisymmetricEquilibrium(g4)

@testset "Bt-Ip Signs" begin
    @test (sign(g1.bcentr), sign(g1.current)) == (1, 1)
    @test (sign(g2.bcentr), sign(g2.current)) == (-1, 1)
    @test (sign(g3.bcentr), sign(g3.current)) == (1, -1)
    @test (sign(g4.bcentr), sign(g4.current)) == (-1, -1)

    @test M1.sigma == 1
    @test M2.sigma == -1
    @test M3.sigma == -1
    @test M4.sigma == 1
end
