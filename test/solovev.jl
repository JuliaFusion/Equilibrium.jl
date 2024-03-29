# ITER parameters
δ = 0.33
ϵ = 0.32
κ = 1.7
B0 = 5.3
R0 = 6.2
qstar = 1.57
alpha = -0.155

S1 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=1,Ip_dir=1)
S2 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=-1,Ip_dir=1)
S3 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=1,Ip_dir=-1)
S4 = solovev(B0, R0, ϵ, δ, κ, alpha, qstar,B0_dir=-1,Ip_dir=-1)

btip1 = (1,1)
btip2 = (-1,1)
btip3 = (1,-1)
btip4 = (-1,-1)

cc1 = 3
cc2 = 3
cc3 = 3
cc4 = 3

r1 = magnetic_axis(S1) .+ (0.1,0.0)
r2 = magnetic_axis(S2) .+ (0.1,0.0)
r3 = magnetic_axis(S3) .+ (0.1,0.0)
r4 = magnetic_axis(S4) .+ (0.1,0.0)

test_data = ((cc1,S1,r1,btip1), (cc2,S2,r2,btip2), (cc3,S3,r3,btip3), (cc4,S4,r4,btip4))

@testset "Solov'ev Tests" begin

    @testset "Bt-Ip Signs" begin
        for (cc, S, r, btip) in test_data
            Bt, Ip = (sign(Bfield(S,r...)[2]), sign(Jfield(S,r...)[2]))
            @test (Bt, Ip) == btip
        end
    end

    @testset verbose = true "CurlB == μ₀J" begin
        for (cc, S, r, btip) in test_data
            @test curlB(S,r...) ≈ mu0*Jfield(S,r...) rtol=0.01
            @test Equilibrium.curlB_autodiff(S,r...) ≈ mu0*Jfield(S,r...) rtol=0.01
        end
    end

    @testset verbose = true "3D tests" begin
        for (cc, S, r, btip) in test_data
            @test Bfield(S,r...) ≈ Bfield(S,r[1],zero(r[1]),r[2]) rtol=0.01
            @test Jfield(S,r...) ≈ Jfield(S,r[1],zero(r[1]),r[2]) rtol=0.01

            @test gradB(S,r...) ≈ gradB(S,r[1],zero(r[1]),r[2]) rtol=0.01
            @test gradB(S,r[1],zero(r[1]),r[2]) ≈ Equilibrium.gradB_autodiff(S,r[1],zero(r[1]),r[2]) rtol=0.01

            @test curlB(S,r...) ≈ curlB(S,r[1],zero(r[1]),r[2]) rtol=0.01
            @test curlB(S,r[1],zero(r[1]),r[2]) ≈ Equilibrium.curlB_autodiff(S,r[1],zero(r[1]),r[2]) rtol=0.01
        end
    end

    @testset verbose = true "COCOS Consistancy" begin
        for (cc, S, r, btip) in test_data
            psi = S(r...)
            sigma_F = sign(poloidal_current(S,psi))
            sigma_pprime = sign(pressure_gradient(S,psi))
            sigma_q = sign(safety_factor(S, psi))
            sigma_dpsi = -sign(psi) # Solov'ev psi at boundary zero
            @test check_cocos(btip[1],btip[2], sigma_F, sigma_pprime, sigma_q, sigma_dpsi, cc; verbose=true)
        end
    end

    @testset "ForwardDiff" begin
        function opti(alpha)
            S = solovev(B0, R0, ϵ, δ, κ, alpha, qstar, B0_dir=1, Ip_dir=1)
            return (S.beta_p - 2)^2
        end

        @test (ForwardDiff.derivative(opti, alpha) !== nothing)
    end

end
