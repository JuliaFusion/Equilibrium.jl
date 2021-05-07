# “One size fits all” analytic solutions to the Grad–Shafranov equation
# Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
#
# Ideal MHD Freidberg Chapter 6.6.1

function _solovev_polynomials(x,y)
    # Eq. 8 & Eq. 27
    ψ₁  = 1
    ψ₂  = x^2
    ψ₃  = y^2 - x^2*log(x)
    ψ₄  = x^4 - 4*x^2*y^2
    ψ₅  = 2*y^4 - 9*y^2*x^2 + 3*x^4*log(x) - 12*x^2*y^2*log(x)
    ψ₆  = x^6 - 12*x^4*y^2 + 8*x^2*y^4
    ψ₇  = 8*y^6 - 140*y^4*x^2 + 75*y^2*x^4 - 15*x^6*log(x) + 180*x^4*y^2*log(x) - 120*x^2*y^4*log(x)
    ψ₈  = y
    ψ₉  = y*x^2
    ψ₁₀ = y^3 - 3*y*x^2*log(x)
    ψ₁₁ = 3*y*x^4 - 4*y^3*x^2
    ψ₁₂ = 8*y^5 - 45*y*x^4 - 80*y^3*x^2*log(x) + 60*y*x^4*log(x)

    ψ_polys = [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

function _solovev_polynomials_Dx(x,y)
    # Derivatives from Mathematica
    ψ₁  = zero(x)
    ψ₂  = 2*x
    ψ₃  = -(2*x*log(x)+ x)
    ψ₄  = 4*x^3 - 8*x*y^2
    ψ₅  = 3*x^3 - 30*x*y^2 + 12*x^3*log(x) - 24*x*y^2*log(x)
    ψ₆  = 6*x^5 - 48*x^3*y^2 + 16*x*y^4
    ψ₇  = -5*x*(3*x^4 - 96*x^2*y^2 + 80*y^4 + 6*(3*x^4 - 24*x^2*y^2 + 8*y^4)*log(x))
    ψ₈  = zero(x)
    ψ₉  = 2*y*x
    ψ₁₀ = -3*x*y*(1 + 2*log(x))
    ψ₁₁ = 12*x^3*y - 8*x*y^3
    ψ₁₂ = 40*x*y*(-3*x^2 - 2*y^2 + (6*x^2 - 4*y^2)*log(x))

    ψ_polys = [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

function _solovev_polynomials_Dxx(x,y)
    # Derivatives from Mathematica
    ψ₁  = zero(x)
    ψ₂  = 2*one(x)
    ψ₃  = -2*log(x) - 3
    ψ₄  = 12*x^2 - 8*y^2
    ψ₅  = 3*(7*x^2 - 18*y^2 + 4*(3*x^2 - 2*y^2)*log(x))
    ψ₆  = 2*(15*x^4 - 72*x^2*y^2 + 8*y^4)
    ψ₇  = -5*(33*x^4 - 432*x^2*y^2 + 128*y^4 + 6*(15*x^4 - 72*x^2*y^2 + 8*y^4)*log(x))
    ψ₈  = zero(x)
    ψ₉  = 2*y
    ψ₁₀ = -3*y*(3 + 2*log(x))
    ψ₁₁ = 36*x^2*y - 8*y^3
    ψ₁₂ = -40*y*(3*x^2 + 6*y^2 - 18*x^2*log(x) + 4*y^2*log(x))

    ψ_polys = [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

function _solovev_polynomials_Dy(x,y)
    # Derivatives from Mathematica
    ψ₁  = zero(y)
    ψ₂  = zero(y)
    ψ₃  = 2*y
    ψ₄  = -8*x^2*y
    ψ₅  = 2*y*(-9*x^2 + 4*y^2 - 12*x^2*log(x))
    ψ₆  = -24*x^4*y + 32*x^2*y^3
    ψ₇  = 2*y*(75*x^4 - 280*x^2*y^2 + 24*y^4 + 60*(3*x^4 - 4*x^2*y^2)*log(x))
    ψ₈  = one(y)
    ψ₉  = x^2
    ψ₁₀ = 3*(y^2 - x^2*log(x))
    ψ₁₁ = 3*(x^4 - 4*x^2*y^2)
    ψ₁₂ = 5*(-9*x^4 + 8*y^4 + 12*(x^4 - 4*x^2*y^2)*log(x))

    ψ_polys = [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

function _solovev_polynomials_Dyy(x,y)
    # Derivatives from Mathematica
    ψ₁  = zero(y)
    ψ₂  = zero(y)
    ψ₃  = 2*one(y)
    ψ₄  = -8*x^2
    ψ₅  = -6*(3*x^2 - 4*y^2 + 4*x^2*log(x))
    ψ₆  = -24*x^2*(x^2 - 4*y^2)
    ψ₇  = 30*(5*x^4 - 56*x^2*y^2 + 8*y^4 + 12*(x^4 - 4*x^2*y^2)*log(x))
    ψ₈  = zero(y)
    ψ₉  = zero(y)
    ψ₁₀ = 6*y
    ψ₁₁ = -24*x^2*y
    ψ₁₂ = 160*y*(y^2 - 3*x^2*log(x))

    ψ_polys = [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

function _solovev_psi_p(x,y,A)
    return x^4/8 + A*(x^2*log(x)/2 - x^4/8)
end

function _solovev_psi_p_Dx(x,y,A)
    # Derivatives from Mathematica
    return (x*(A + x^2 - A*x^2 + 2*A*log(x)))/2
end

function _solovev_psi_p_Dxx(x,y,A)
    # Derivatives from Mathematica
    return (3*(A + x^2 - A*x^2))/2 + A*log(x)
end

function _solovev_psi_p_Dy(x,y,A)
    # Derivatives from Mathematica
    return zero(y)
end

function _solovev_psi_p_Dyy(x,y,A)
    # Derivatives from Mathematica
    return zero(y)
end

function _solovev_psi(x,y,A,coeffs)

    ψ_p = _solovev_psi_p(x,y,A)
    ψ_polys = _solovev_polynomials(x,y)

    ψ_h = zero(x)
    for i in eachindex(coeffs)
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dx(x,y,A,coeffs)

    ψ_p = _solovev_psi_p_Dx(x,y,A)
    ψ_polys = _solovev_polynomials_Dx(x,y)

    ψ_h = zero(x)
    for i in eachindex(coeffs)
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dxx(x,y,A,coeffs)

    ψ_p = _solovev_psi_p_Dxx(x,y,A)
    ψ_polys = _solovev_polynomials_Dxx(x,y)

    ψ_h = zero(x)
    for i in eachindex(coeffs)
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dy(x,y,A,coeffs)

    ψ_p = _solovev_psi_p_Dy(x,y,A)
    ψ_polys = _solovev_polynomials_Dy(x,y)

    ψ_h = zero(y)
    for i in eachindex(coeffs)
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dyy(x,y,A,coeffs)

    ψ_p = _solovev_psi_p_Dyy(x,y,A)
    ψ_polys = _solovev_polynomials_Dyy(x,y)

    ψ_h = zero(x)
    for i in eachindex(coeffs)
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

"""
SolovevEquilibrium Structure

Defines the equilibrium that satisfy Δ⋆ψ(x,y) = α + (1 - α)x^2 where r,z = R0*(x,y)

`cocos` - COCOS

`R0` - Major Radius [m]

`epsilon` - inverse aspect ratio a/R0

`delta` - triangularity

`kappa` - elongation/ellipticity

`alpha` - constant relating beta regime

`qstar` - Kink safety factor

`psi0` - Poloidal flux normalization

`beta_p` - Poloidal beta

`beta_t` - Toroidal beta

`xp` - x-point

`c` - Coefficients for Solov'ev polynomials

`symmetric` - If true then equilibrium is up-down symmetric

"""
struct SolovevEquilibrium
    cocos::COCOS
    B0
    R0
    epsilon
    delta
    kappa
    alpha
    qstar
    psi0
    beta_p
    beta_t
    xp::Union{NTuple{2},Nothing}
    c::AbstractVector
    symmetric::Bool
end

function SolovevEquilibrium(B0, R0, ϵ, δ, κ, α, q⋆;
                            diverted::Bool = false,
                            xpoint::Union{NTuple{2,<:Real}} = diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing,
                            symmetric::Bool = (xpoint == nothing))

    if δ > sin(1)
        @warn "Equilibrium is not convex. δ > sin(1)"
    end

    #Eq. 11
    δ₀ = asin(δ)
    N₁ = -(1 + δ₀)^2/(ϵ*κ^2)
    N₂ =  (1 - δ₀)^2/(ϵ*κ^2)
    N₃ = -κ/(ϵ*cos(δ₀)^2)

    Ψ_h = zeros(12,12)
    Ψ_p = zeros(12)
    # Boundary Conditions: (ψ_p - bc) + ψ_h'*c = 0
    if symmetric
        if xpoint == nothing
            # Eq. 10
            # Outer equatorial point
            x, y = 1 + ϵ, 0
            Ψ_p[1] = _solovev_psi_p(x, y, α)
            Ψ_h[:,1] .= _solovev_polynomials(x,y)
            # Inner equatorial point
            x, y = 1 - ϵ, 0
            Ψ_p[2] = _solovev_psi_p(x, y, α)
            Ψ_h[:,2] .= _solovev_polynomials(x, y)
            # High point
            x, y = 1 - δ*ϵ, κ*ϵ
            Ψ_p[3] = _solovev_psi_p(x, y, α)
            Ψ_h[:,3] .= _solovev_polynomials(x, y)
            # High point maximum
            x, y = 1 - δ*ϵ, κ*ϵ
            Ψ_p[4] = _solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,4] .= _solovev_polynomials_Dx(x, y)
            # Outer equatorial point curvature
            x, y = 1 + ϵ, 0
            Ψ_p[5] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_p_Dx(x, y, α)
            Ψ_h[:,5] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
            # Inner equatorial point curvature
            x, y = 1 - ϵ, 0
            Ψ_p[6] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_p_Dx(x, y, α)
            Ψ_h[:,6] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)
            # High point curvature
            x, y = 1 - δ*ϵ, κ*ϵ
            Ψ_p[7] = _solovev_psi_p_Dxx(x, y, α) + N₃*_solovev_p_Dy(x, y, α)
            Ψ_h[:,7] .= _solovev_polynomials_Dxx(x, y) .+ N₃*_solovev_polynomials_Dy(x, y)

            c = Ψ_h[1:7,1:7]'\(-Ψ_p[1:7])
        else
            xsep, ysep = xpoint[1]/R0, xpoint[2]/R0 # normalize xpoint

            # Eq. 12
            # Outer equatorial point
            x, y = 1 + ϵ, 0
            Ψ_p[1] = _solovev_psi_p(x, y, α)
            Ψ_h[:,1] .= _solovev_polynomials(x,y)
            # Inner equatorial point
            x, y = 1 - ϵ, 0
            Ψ_p[2] = _solovev_psi_p(x, y, α)
            Ψ_h[:,2] .= _solovev_polynomials(x, y)
            # High point
            x, y = xsep, ysep
            Ψ_p[3] = _solovev_psi_p(x, y, α)
            Ψ_h[:,3] .= _solovev_polynomials(x, y)
            # B_norm = 0 at high point
            x, y = xsep, ysep
            Ψ_p[4] = _solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,4] .= _solovev_polynomials_Dx(x, y)
            # B_tang = 0 at high point
            x, y = xsep, ysep
            Ψ_p[5] = _solovev_psi_p_Dy(x, y, α)
            Ψ_h[:,5] .= _solovev_polynomials_Dy(x, y)
            # Outer equatorial point curvature
            x, y = 1 + ϵ, 0
            Ψ_p[6] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_p_Dx(x, y, α)
            Ψ_h[:,6] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
            # Inner equatorial point curvature
            x, y = 1 - ϵ, 0
            Ψ_p[7] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_p_Dx(x, y, α)
            Ψ_h[:,7] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)

            c = Ψ_h[1:7,1:7]'\(-Ψ_p[1:7])
        end
    else
        xsep, ysep = xpoint[1]/R0, xpoint[2]/R0 # normalize xpoint
        if ysep > 0
            throw(ArgumentError("X-point should be below the midplane"))
        end

        # Eq. 28
        # Outer equatorial point
        x, y = 1 + ϵ, 0
        Ψ_p[1] = _solovev_psi_p(x, y, α)
        Ψ_h[:,1] .= _solovev_polynomials(x,y)
        # Inner equatorial point
        x, y = 1 - ϵ, 0
        Ψ_p[2] = _solovev_psi_p(x, y, α)
        Ψ_h[:,2] .= _solovev_polynomials(x, y)
        # Upper high point
        x, y = 1 - δ*ϵ, κ*ϵ
        Ψ_p[3] = _solovev_psi_p(x, y, α)
        Ψ_h[:,3] .= _solovev_polynomials(x, y)
        # Lower X-point
        x, y = xsep, ysep
        Ψ_p[4] = _solovev_psi_p(x, y, α)
        Ψ_h[:,4] .= _solovev_polynomials(x, y)
        # Outer equatorial point up-down symmetry
        x, y = 1 + ϵ, 0
        Ψ_p[5] = _solovev_psi_p_Dy(x,y, α)
        Ψ_h[:,5] .= _solovev_polynomials_Dy(x,y)
        # Inner equatorial point up-down symmetry
        x, y = 1 - ϵ, 0
        Ψ_p[6] = _solovev_psi_p_Dy(x, y, α)
        Ψ_h[:,6] .= _solovev_polynomials_Dy(x, y)
        # Upper high point maximum
        x, y = 1 - δ*ϵ, κ*ϵ
        Ψ_p[7] = _solovev_psi_p_Dx(x, y, α)
        Ψ_h[:,7] .= _solovev_polynomials_Dx(x, y)
        # B_y = 0 at X-point
        x, y = xsep, ysep
        Ψ_p[8] = _solovev_psi_p_Dx(x, y, α)
        Ψ_h[:,8] .= _solovev_polynomials_Dx(x, y)
        # B_x = 0 at X-point
        x, y = xsep, ysep
        Ψ_p[9] = _solovev_psi_p_Dy(x, y, α)
        Ψ_h[:,9] .= _solovev_polynomials_Dy(x, y)
        # Outer equatorial point curvature
        x, y = 1 + ϵ, 0
        Ψ_p[10] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_p_Dx(x, y, α)
        Ψ_h[:,10] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
        # Inner equatorial point curvature
        x, y = 1 - ϵ, 0
        Ψ_p[11] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_p_Dx(x, y, α)
        Ψ_h[:,11] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)
        # High point curvature
        x, y = 1 - δ*ϵ, κ*ϵ
        Ψ_p[12] = _solovev_psi_p_Dxx(x, y, α) + N₃*_solovev_p_Dy(x, y, α)
        Ψ_h[:,12] .= _solovev_polynomials_Dxx(x, y) .+ N₃*_solovev_polynomials_Dy(x, y)

        c = Ψ_h'\(-Ψ_p)
    end

    τ = range(0,2pi,length=100)
    x = (1 .+ ϵ .* cos.(τ .+ δ₀*sin.(τ)))
    y = (ϵ*κ*sin.(τ))
    bdry = PlasmaBoundary(collect(zip(x,y)))

    #Eq. 6.158 in Ideal MHD
    K₁ = area_average(bdry, (x,y) -> (α + (1-α)*x^2)/x)
    K₂ = volume_average(bdry,(x,y) -> -_solovev_psi(x,y,α,c))/(2pi)
    K₃ = volume(bdry)/(2pi)

    #Eq. 6.157 in Ideal MHD
    a = ϵ*R0
    ψ0 = 2pi*(a^2*B0)*(1 + κ^2)/(2*q⋆*K₁)
    β_p = 8*pi^2*ϵ^2*(1-α)*((1 + κ^2)/2)*K₂/(K₁^2*K₃)
    β_t = (8*pi^2*ϵ^4*(1-α)/(q⋆^2))*(((1+κ^2)/2)^2)*K₂/(K₁^2*K₃)

    return SolovevEquilibrium(cocos(3),B0,R0,ϵ,δ,κ,α,q⋆,ψ0,β_p,β_t,xpoint,c,symmetric)
end

function (S::SolovevEquilibrium)(r,z)
    x = r/S.R0
    y = z/S.R0
    return S.psi0*_solovev_psi(x,y,S.alpha,S.c)
end

function Base.show(io::IO, S::SolovevEquilibrium)
    print(io, "SolovevEquilibrium\n")
    print(io, "  B0 = $(S.B0) [T]\n")
    print(io, "  R0 = $(S.R0) [m]\n")
    print(io, "  ϵ  = $(S.epsilon)\n")
    print(io, "  δ  = $(S.delta)\n")
    print(io, "  κ  = $(S.kappa)\n")
    print(io, "  α  = $(S.alpha)\n")
    print(io, "  q⋆ = $(S.qstar)\n")
    print(io, "  βp = $(S.beta_p)\n")
    print(io, "  βt = $(S.beta_t)\n")
end

Base.broadcastable(S::SolovevEquilibrium) = (S,)
