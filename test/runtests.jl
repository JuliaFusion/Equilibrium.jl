using Equilibrium
using Test
using EFIT

dir = @__DIR__
g1 = readg(dir*"/g150219.03200") # (+Bt, +Ip)
g2 = readg(dir*"/g133221.01150") # (-Bt, +Ip)
g3 = readg(dir*"/g159177.02700") # (+Bt, -Ip)
g4 = readg(dir*"/g153298.04400") # (-Bt, -Ip)

# COCOS for the GEQDSK files
cc1 = 1
cc2 = 5
cc3 = 3
cc4 = 7

M1 = AxisymmetricEquilibrium(g1)
M2 = AxisymmetricEquilibrium(g2)
M3 = AxisymmetricEquilibrium(g3)
M4 = AxisymmetricEquilibrium(g4)

r1 = (M1.axis[1]+0.1,M1.axis[2])
r2 = (M2.axis[1]+0.1,M2.axis[2])
r3 = (M3.axis[1]+0.1,M3.axis[2])
r4 = (M4.axis[1]+0.1,M4.axis[2])

@testset verbose = true "Equilibrium Tests" begin
    @testset "Bt-Ip Signs" begin
        @test (sign(g1.bcentr), sign(g1.current)) == (1, 1)
        @test (sign(g2.bcentr), sign(g2.current)) == (-1, 1)
        @test (sign(g3.bcentr), sign(g3.current)) == (1, -1)
        @test (sign(g4.bcentr), sign(g4.current)) == (-1, -1)

        Bt1, Ip1 = (sign(Bfield(M1,r1...)[2]), sign(Jfield(M1,r1...)[2]))
        @test (Bt1, Ip1) == (1, 1)

        (Bt2, Ip2) = (sign(Bfield(M2,r2...)[2]), sign(Jfield(M2,r2...)[2]))
        @test (Bt2, Ip2) == (-1, 1)

        (Bt3, Ip3) = (sign(Bfield(M3,r3...)[2]), sign(Jfield(M3,r2...)[2]))
        @test (Bt3, Ip3) == (1, -1)

        (Bt4, Ip4) = (sign(Bfield(M4,r4...)[2]), sign(Jfield(M4,r4...)[2]))
        @test (Bt4, Ip4) == (-1, -1)
    end

    @testset verbose = true "CurlB == μ₀J" begin
        mu0 = 4pi*10^-7
        @test curlB(M1,r1...) ≈ mu0*Jfield(M1,r1...) rtol=0.01
        @test curlB(M2,r2...) ≈ mu0*Jfield(M2,r2...) rtol=0.01
        @test curlB(M3,r3...) ≈ mu0*Jfield(M3,r3...) rtol=0.01
        @test curlB(M4,r4...) ≈ mu0*Jfield(M4,r4...) rtol=0.01
    end

    @testset verbose = true "COCOS Identify" begin
        @test identify_cocos(g1; clockwise_phi = false) == (cc1,)
        @test identify_cocos(g2; clockwise_phi = false) == (cc2,)
        @test identify_cocos(g3; clockwise_phi = false) == (cc3,)
        @test identify_cocos(g4; clockwise_phi = false) == (cc4,)
    end

    @testset verbose = true "COCOS Consistancy" begin
        @test check_cocos(g1, cc1)
        @test check_cocos(g2, cc2)
        @test check_cocos(g3, cc3)
        @test check_cocos(g4, cc4)
    end

    @testset verbose = true "COCOS GEQDSK Transformation" begin
        for (cc_in, g) in ((cc1, g1), (cc2,g2), (cc3, g3), (cc4, g4))
            for cc = 1:8
                _cc = cocos(cc)
                g_new = transform_cocos(g, cc_in, _cc)
                @test identify_cocos(g_new, clockwise_phi = _cc.sigma_RpZ < 0) == (cc,)
                @test check_cocos(g_new, cc)
            end
        end
    end
end
