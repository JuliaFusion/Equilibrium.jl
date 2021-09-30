# “One size fits all” analytic solutions to the Grad–Shafranov equation
# Phys. Plasmas 17, 032502 (2010); https://doi.org/10.1063/1.3328818
#
# Ideal MHD Freidberg Chapter 6.6.1

@fastmath function _solovev_polynomials(x, y, ::Val{N} = Val(12); logx=log(x)) where N

    x2 = x^2
    x4 = x2^2
    x6 = x2^3
    y2 = y^2
    y4 = y2^2
    y6 = y2^3

    # Eq. 8 & Eq. 27
    ψ₁  = 1
    ψ₂  = x2
    ψ₃  = y2 - x2*logx
    ψ₄  = x4 - 4*x2*y2
    ψ₅  = 2*y4 - 9*y2*x2 + 3*x4*logx - 12*x2*y2*logx
    ψ₆  = x6 - 12*x4*y2 + 8*x2*y4
    ψ₇  = 8*y6 - 140*y4*x2 + 75*y2*x4 - 15*x6*logx + 180*x4*y2*logx - 120*x2*y4*logx

    N == 7 && return @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇]

    x3 = x^3
    y3 = y^3
    y5 = y^5

    ψ₈  = y
    ψ₉  = y*x2
    ψ₁₀ = y3 - 3*y*x2*logx
    ψ₁₁ = 3*y*x4 - 4*y3*x2
    ψ₁₂ = 8*y5 - 45*y*x4 - 80*y3*x2*logx + 60*y*x4*logx

    ψ_polys = @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

@fastmath function _solovev_polynomials_Dx(x, y, ::Val{N} = Val(12); logx=log(x)) where N

    x2 = x^2
    x3 = x^3
    x4 = x2^2
    x5 = x^5
    x6 = x2^3
    y2 = y^2
    y4 = y2^2

    # Derivatives from Mathematica
    ψ₁  = zero(x)
    ψ₂  = 2*x
    ψ₃  = -(2*x*logx+ x)
    ψ₄  = 4*x3 - 8*x*y2
    ψ₅  = 3*x3 - 30*x*y2 + 12*x3*logx - 24*x*y2*logx
    ψ₆  = 6*x5 - 48*x3*y2 + 16*x*y4
    ψ₇  = -5*x*(3*x4 - 96*x2*y2 + 80*y4 + 6*(3*x4 - 24*x2*y2 + 8*y4)*logx)

    N == 7 && return @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇]

    y3 = y^3

    ψ₈  = zero(x)
    ψ₉  = 2*y*x
    ψ₁₀ = -3*x*y*(1 + 2*logx)
    ψ₁₁ = 12*x3*y - 8*x*y3
    ψ₁₂ = 40*x*y*(-3*x2 - 2*y2 + (6*x2 - 4*y2)*logx)

    ψ_polys = @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

@fastmath function _solovev_polynomials_Dxx(x, y, ::Val{N} = Val(12); logx=log(x)) where N

    x2 = x^2
    x3 = x^3
    x4 = x2^2
    x6 = x2^3
    y2 = y^2
    y4 = y2^2

    # Derivatives from Mathematica
    ψ₁  = zero(x)
    ψ₂  = 2*one(x)
    ψ₃  = -2*logx - 3
    ψ₄  = 12*x2 - 8*y2
    ψ₅  = 3*(7*x2 - 18*y2 + 4*(3*x2 - 2*y2)*logx)
    ψ₆  = 2*(15*x4 - 72*x2*y2 + 8*y4)
    ψ₇  = -5*(33*x4 - 432*x2*y2 + 128*y4 + 6*(15*x4 - 72*x2*y2 + 8*y4)*logx)

    N == 7 && return @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇]

    y3 = y^3

    ψ₈  = zero(x)
    ψ₉  = 2*y
    ψ₁₀ = -3*y*(3 + 2*logx)
    ψ₁₁ = 36*x2*y - 8*y3
    ψ₁₂ = -40*y*(3*x2 + 6*y2 - 18*x2*logx + 4*y2*logx)

    ψ_polys = @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

@fastmath function _solovev_polynomials_Dy(x, y, ::Val{N} = Val(12); logx=log(x)) where N

    x2 = x^2
    x3 = x^3
    x4 = x2^2
    x6 = x2^3
    y2 = y^2
    y3 = y^3
    y4 = y2^2

    # Derivatives from Mathematica
    ψ₁  = zero(y)
    ψ₂  = zero(y)
    ψ₃  = 2*y
    ψ₄  = -8*x2*y
    ψ₅  = 2*y*(-9*x2 + 4*y2 - 12*x2*logx)
    ψ₆  = -24*x4*y + 32*x2*y3
    ψ₇  = 2*y*(75*x4 - 280*x2*y2 + 24*y4 + 60*(3*x4 - 4*x2*y2)*logx)

    N == 7 && return @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇]

    ψ₈  = one(y)
    ψ₉  = x2
    ψ₁₀ = 3*(y2 - x2*logx)
    ψ₁₁ = 3*(x4 - 4*x2*y2)
    ψ₁₂ = 5*(-9*x4 + 8*y4 + 12*(x4 - 4*x2*y2)*logx)

    ψ_polys = @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

@fastmath function _solovev_polynomials_Dyy(x, y, ::Val{N} = Val(12); logx=log(x)) where N

    x2 = x^2
    x3 = x^3
    x4 = x2^2
    x6 = x2^3
    y2 = y^2
    y3 = y^3
    y4 = y2^2

    # Derivatives from Mathematica
    ψ₁  = zero(y)
    ψ₂  = zero(y)
    ψ₃  = 2*one(y)
    ψ₄  = -8*x2
    ψ₅  = -6*(3*x2 - 4*y2 + 4*x2*logx)
    ψ₆  = -24*x2*(x2 - 4*y2)
    ψ₇  = 30*(5*x4 - 56*x2*y2 + 8*y4 + 12*(x4 - 4*x2*y2)*logx)

    N == 7 && return @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇]

    ψ₈  = zero(y)
    ψ₉  = zero(y)
    ψ₁₀ = 6*y
    ψ₁₁ = -24*x2*y
    ψ₁₂ = 160*y*(y2 - 3*x2*logx)

    ψ_polys = @SVector [ψ₁, ψ₂, ψ₃, ψ₄, ψ₅, ψ₆, ψ₇, ψ₈, ψ₉, ψ₁₀, ψ₁₁, ψ₁₂]

    return ψ_polys
end

@fastmath function _solovev_psi_p(x,y,A; logx=log(x))
    # Derivatives from Mathematica
    x2 = x^2
    x4 = x2*x2
    return x4/8 + A*(x2*logx/2 - x4/8)
end

@fastmath function _solovev_psi_p_Dx(x,y,A; logx=log(x))
    # Derivatives from Mathematica
    x2 = x^2
    return (x*(A + x2 - A*x2 + 2*A*logx))/2
end

@fastmath function _solovev_psi_p_Dxx(x,y,A; logx=log(x))
    # Derivatives from Mathematica
    x2 = x^2
    return (3*(A + x2 - A*x2))/2 + A*logx
end

function _solovev_psi_p_Dy(x,y,A)
    return zero(y)
end

function _solovev_psi_p_Dyy(x,y,A)
    return zero(y)
end

function _solovev_psi(x,y,A,coeffs::SVector{N};logx=log(x)) where N

    ψ_p = _solovev_psi_p(x,y,A; logx=logx)
    ψ_polys = _solovev_polynomials(x,y,Val(N),logx=logx)

    ψ_h = zero(x)
    @inbounds for i=1:N
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dx(x,y,A,coeffs::SVector{N};logx=log(x)) where N

    ψ_p = _solovev_psi_p_Dx(x,y,A,logx=logx)
    ψ_polys = _solovev_polynomials_Dx(x,y,Val(N),logx=logx)

    ψ_h = zero(x)
    @inbounds for i=1:N
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dxx(x,y,A,coeffs::SVector{N};logx=log(x)) where N

    ψ_p = _solovev_psi_p_Dxx(x,y,A,logx=logx)
    ψ_polys = _solovev_polynomials_Dxx(x,y,Val(N),logx=logx)

    ψ_h = zero(x)
    @inbounds for i=1:N
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dy(x,y,A,coeffs::SVector{N};logx=log(x)) where N

    ψ_p = _solovev_psi_p_Dy(x,y,A)
    ψ_polys = _solovev_polynomials_Dy(x,y,Val(N),logx=logx)

    ψ_h = zero(y)
    @inbounds for i=1:N
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

function _solovev_psi_Dyy(x,y,A,coeffs::SVector{N};logx=log(x)) where N

    ψ_p = _solovev_psi_p_Dyy(x,y,A)
    ψ_polys = _solovev_polynomials_Dyy(x,y,Val(N),logx=logx)

    ψ_h = zero(x)
    @inbounds for i=1:N
        ψ_h += coeffs[i]*ψ_polys[i]
    end

    return ψ_p + ψ_h
end

"""
SolovevEquilibrium Structure

Defines the equilibrium that satisfy Δ⋆ψ(x,y) = α + (1 - α)x^2\\
where (r,z) = (R0 x,R0 y). F(dF/dψ) = -A, μ₀dp/dψ = -C, α = A/(A + C R0²), ψ₀ = R0²(A + C R0²)

Fields:\\
`cocos` - COCOS\\
`B0` - Toroidal magnetic field magnitude on axis [T]\\
`R0` - Major Radius [m]\\
`epsilon` - inverse aspect ratio a/R0\\
`delta` - triangularity\\
`kappa` - elongation/ellipticity\\
`alpha` - constant relating beta regime (α)\\
`qstar` - Kink safety factor\\
`psi0` - Poloidal flux normalization (ψ₀)\\
`beta_p` - Poloidal beta\\
`beta_t` - Toroidal beta\\
`c` - Coefficients for Solov'ev polynomials\\
`diverted` - If true then equilibrium has one or more x-points\\
`symmetric` - If true then equilibrium is up-down symmetric\\
`sigma` - Sign of dot(B,J)
"""
struct SolovevEquilibrium{T,N} <: AbstractEquilibrium
    cocos::COCOS
    B0::T
    R0::T
    epsilon::T
    delta::T
    kappa::T
    alpha::T
    qstar::T
    psi0::T
    beta_p::T
    beta_t::T
    c::SVector{N,T}
    diverted::Bool
    symmetric::Bool
    sigma_B0::Int
    sigma_Ip::Int
end

"""
    solovev(B0, R0, ϵ, δ, κ, q⋆; B0_dir=1, Ip_dir=1, diverted=false, xpoint=nothing, symmetric = true)

Creates a SolovevEquilibrium Structure

Arguments:\\
`B0` - Toroidal magnetic field magnitude on axis [T]\\
`R0` - Major Radius [m]\\
`ϵ` - inverse aspect ratio a/R0\\
`δ` - triangularity\\
`κ` - elongation/ellipticity\\
`α` - constant relating beta regime (α)\\
`q⋆` - Kink safety factor\\

Keyword Arguments:\\
`B0_dir` - Direction of the toroidal field (+1 = CCW, -1 = CW)\\
`Ip_dir` - Direction of the plasma current (+1 = CCW, -1 = CW)\\
`diverted` - If true then equilibrium has one or more x-points\\
`xpoint` - If diverted = true, then xpoint is set to (R0*(1-1.1*δ*ϵ), -R0*1.1*κ*ϵ) else nothing\\
`symmetric` - Is equilibrium up-down symmetric\\
"""
function solovev(B0, R0, ϵ, δ, κ, α, qstar;
                 B0_dir = 1, Ip_dir = 1,
                 diverted::Bool = false,
                 xpoint::Union{NTuple{2},Nothing} = (diverted ? (R0*(1-1.1*δ*ϵ),-R0*1.1*κ*ϵ) : nothing),
                 symmetric::Bool = (xpoint === nothing))

    if δ > sin(1)
        @warn "Equilibrium is not convex. δ > sin(1)"
    end

    diverted = xpoint !== nothing

    #Eq. 11
    δ₀ = asin(δ)
    N₁ = -(1 + δ₀)^2/(ϵ*κ^2) # [d²x/dy²]_(τ=0)
    N₂ =  (1 - δ₀)^2/(ϵ*κ^2) # [d²x/dy²]_(τ=pi)
    N₃ = -κ/(ϵ*cos(δ₀)^2)    # [d²y/dx²]_(τ=pi/2)

    Ψ_h = zeros(typeof(α),12,12)
    Ψ_p = zeros(typeof(α),12)
    # Boundary Conditions: (ψ_p - bc) + ψ_h'*c = 0
    if symmetric
        if xpoint === nothing
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
            Ψ_p[5] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,5] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
            # Inner equatorial point curvature
            x, y = 1 - ϵ, 0
            Ψ_p[6] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,6] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)
            # High point curvature
            x, y = 1 - δ*ϵ, κ*ϵ
            Ψ_p[7] = _solovev_psi_p_Dxx(x, y, α) + N₃*_solovev_psi_p_Dy(x, y, α)
            Ψ_h[:,7] .= _solovev_polynomials_Dxx(x, y) .+ N₃*_solovev_polynomials_Dy(x, y)

            c = SVector{7}(Ψ_h[1:7,1:7]'\(-Ψ_p[1:7]))
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
            Ψ_p[6] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,6] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
            # Inner equatorial point curvature
            x, y = 1 - ϵ, 0
            Ψ_p[7] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_psi_p_Dx(x, y, α)
            Ψ_h[:,7] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)

            c = SVector{7}(Ψ_h[1:7,1:7]'\(-Ψ_p[1:7]))
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
        Ψ_p[10] = _solovev_psi_p_Dyy(x, y, α) + N₁*_solovev_psi_p_Dx(x, y, α)
        Ψ_h[:,10] .= _solovev_polynomials_Dyy(x, y) .+ N₁*_solovev_polynomials_Dx(x, y)
        # Inner equatorial point curvature
        x, y = 1 - ϵ, 0
        Ψ_p[11] = _solovev_psi_p_Dyy(x, y, α) + N₂*_solovev_psi_p_Dx(x, y, α)
        Ψ_h[:,11] .= _solovev_polynomials_Dyy(x, y) .+ N₂*_solovev_polynomials_Dx(x, y)
        # High point curvature
        x, y = 1 - δ*ϵ, κ*ϵ
        Ψ_p[12] = _solovev_psi_p_Dxx(x, y, α) + N₃*_solovev_psi_p_Dy(x, y, α)
        Ψ_h[:,12] .= _solovev_polynomials_Dxx(x, y) .+ N₃*_solovev_polynomials_Dy(x, y)

        c = SVector{12}(Ψ_h'\(-Ψ_p))
    end

    τ = range(0,2pi,length=100)
    x = (1 .+ ϵ .* cos.(τ .+ δ₀*sin.(τ)))
    y = (ϵ*κ*sin.(τ))
    bdry = PlasmaBoundary(collect(zip(x,y)))

    #Eq. 6.158 in Ideal MHD
    V = volume(bdry,dx=0.005,dy=0.005)
    K₁ = area_average(bdry, (x,y) -> (α + (1-α)*x^2)/x,dx=0.005,dy=0.005)*area(bdry)
    K₂ = volume_average(bdry,(x,y) -> -_solovev_psi(x,y,α,c),dx=0.005,dy=0.005)*V/(2pi)
    K₃ = V/(2pi)

    #Eq. 19 from paper (error in β_t in Ideal MHD (extra ^2))
    Cp = circumference(bdry)
    a = ϵ*R0
    ψ0 = -Cp*a*R0*B0/(qstar*K₁*Ip_dir)
    β_p = 2*(1-α)*Cp^2 * K₂/(K₁^2*K₃)
    β_t = ϵ^2*β_p/(qstar^2)

    cc = 3
    return SolovevEquilibrium(cocos(cc),promote(B0,R0,ϵ,δ,κ,α,qstar,ψ0,β_p,β_t)...,c, diverted, symmetric, B0_dir, Ip_dir)
end

function (S::SolovevEquilibrium)(r,z)
    x = r/S.R0
    y = z/S.R0
    logx = log(x)
    return S.psi0*_solovev_psi(x,y,S.alpha,S.c,logx=logx)
end

function Base.show(io::IO, S::SolovevEquilibrium)
    print(io, "SolovevEquilibrium\n")
    print(io, "  B0 = $(S.sigma_B0*S.B0) [T]\n")
    print(io, "  R0 = $(S.R0) [m]\n")
    print(io, "  ϵ  = $(S.epsilon)\n")
    print(io, "  δ  = $(S.delta)\n")
    print(io, "  κ  = $(S.kappa)\n")
    print(io, "  α  = $(S.alpha)\n")
    print(io, "  q⋆ = $(S.qstar)\n")
    print(io, "  βp = $(S.beta_p)\n")
    print(io, "  βt = $(S.beta_t)\n")
    print(io, "  σ  = $(S.sigma_B0*S.sigma_Ip)\n")
    print(io, "  diverted  = $(S.diverted)\n")
    print(io, "  symmetric = $(S.symmetric)")
end

Base.broadcastable(S::SolovevEquilibrium) = (S,)

function cocos(S::SolovevEquilibrium)
    return S.cocos
end

function B0Ip_sign(S::SolovevEquilibrium)
    return S.sigma_B0*S.sigma_Ip
end

function limits(S::SolovevEquilibrium)
    xlims = (S.R0*(1-1.2*S.epsilon), S.R0*(1 + 1.2*S.epsilon))
    ylims = (S.R0*(-1.2*S.epsilon*S.kappa), S.R0*(1.2*S.epsilon*S.kappa))
    return xlims, ylims
end

function boundary(S::SolovevEquilibrium; n=100, kwargs...)
    if S.diverted
        return boundary(S,0.0; kwargs...)
    else
        τ = range(0,2pi,length=n)
        δ₀ = asin(S.delta)
        x = S.R0*(1 .+ S.epsilon .* cos.(τ .+ δ₀*sin.(τ)))
        y = S.R0*(S.epsilon*S.kappa*sin.(τ))
        return PlasmaBoundary(collect(zip(x,y)))
    end
end

function psi_gradient(S::SolovevEquilibrium,r,z)
    x, y = r/S.R0, z/S.R0
    logx=log(x)
    return (S.psi0/S.R0)*SVector{2}(_solovev_psi_Dx(x,y,S.alpha,S.c; logx=logx), _solovev_psi_Dy(x,y,S.alpha,S.c;logx=logx))
end

_solovev_magnetic_axis = Dict{SolovevEquilibrium,NTuple{2}}()
function clear_cache(S)
    delete!(_solovev_magnetic_axis,S)
end

function magnetic_axis(S::SolovevEquilibrium; x0 = (S.R0, zero(S.R0)))
    if S in keys(_solovev_magnetic_axis)
        return _solovev_magnetic_axis[S]
    else
        sigma_psi = sign(S(x0[1],x0[2]))
        res = Optim.optimize(x->-sigma_psi*S(x[1],x[2]), x -> -sigma_psi*psi_gradient(S,x[1],x[2]), collect(x0), inplace=false)
        axis = (res.minimizer[1],res.minimizer[2])
        _solovev_magnetic_axis[S] = axis
    end
    return axis
end

function psi_limits(S::SolovevEquilibrium)
    psimag = S(magnetic_axis(S)...)
    psibry = zero(psimag)
    return (psimag, psibry)
end

function pressure_gradient(S::SolovevEquilibrium)
    C = S.psi0*(1-S.alpha)/(S.R0^4)
    return -C/mu0
end

function pressure_gradient(S::SolovevEquilibrium, psi)
    return pressure_gradient(S) + psi*0
end

function pressure(S::SolovevEquilibrium, psi; p0=zero(psi))
    C = S.psi0*(1-S.alpha)/(S.R0^4)
    dpdpsi = -C/mu0
    p = dpdpsi*psi + p0
    return p
end

function poloidal_current(S::SolovevEquilibrium,psi)
    A = S.psi0*S.alpha/(S.R0^2)
    dF2dpsi = -2*A
    F2 = abs(dF2dpsi*psi + S.R0^2*S.B0^2)
    return S.sigma_B0*sqrt(F2)
end

function poloidal_current_gradient(S::SolovevEquilibrium,psi)
    A = S.psi0*S.alpha/(S.R0^2)
    return -A/poloidal_current(S,psi)
end

electric_potential(S::SolovevEquilibrium,psi) = zero(psi)
electric_potential_gradient(S::SolovevEquilibrium,psi) = zero(psi)

# === Special Cases ===

function curlB(S::SolovevEquilibrium, r, z)
    # curlB = mu0*J
    psi = S(r,z)
    gval = poloidal_current(S,psi)
    grad_psi = psi_gradient(S,r,z)

    gp = poloidal_current_gradient(S,psi)
    pp = pressure_gradient(S)

    cc = cocos(S)

    invr = inv(r)
    cocos_factor = cc.sigma_RpZ*invr
    cBr = -cocos_factor*gp*grad_psi[2]
    cBz =  cocos_factor*gp*grad_psi[1]

    cocos_factor = -cc.sigma_Bp*((2pi)^cc.exp_Bp)
    cBt = cocos_factor*(r*pp*mu0 + gval*gp*invr)

    return SVector{3}(cylindrical_cocos(cc,cBr,cBt,cBz))
end

function curlB(S::SolovevEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    cB = curlB(S,r,z)
    sp, cp = sincos(phi)
    cB_xyz = SVector{3}(cB[1]*cp - cB[2]*sp, cB[1]*sp + cB[2]*cp, cB[3])
    return cB_xyz
end
