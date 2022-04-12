abstract type PlasmaShape end

aspect_ratio(S::PlasmaShape) = inv(S.ϵ)
elongation(S::PlasmaShape) = S.κ
major_radius(S::PlasmaShape) = S.R0
minor_radius(S::PlasmaShape) = S.R0*S.ϵ

function scale_aspect(S::PlasmaShape,s)
    SS = copy(S)
    SS.ϵ = s*S.ϵ
end

"""
MillerShape Structure

Defines the Miller Plasma Shape Parameterization

Fields:\\
`R0` - Major Radius [m]\\
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`δ`  - Triangularity
"""
mutable struct MillerShape{T<:Number} <: PlasmaShape
    R0::T  # Major Radius [m]
    ϵ::T   # Inverse Aspect Ratio a/R0 (a = minor radius)
    κ::T   # Elongation
    δ::T   # Triangularity
end

function MillerShape(R0,ϵ,κ,δ)
    MillerShape(promote(R0,ϵ,κ,δ)...)
end

elevation(S::MillerShape) = zero(S.δ)
triangularity(S::MillerShape) = S.δ
tilt(S::MillerShape) = zero(S.δ)
ovality(S::MillerShape) = zero(S.δ)
squareness(S::MillerShape) = zero(S.δ)

function shape(S::MillerShape; N=100)
    τ = range(0,2pi,length=N)
    δ₀ = asin(S.δ)
    x = S.R0*(1 .+ S.ϵ .* cos.(τ .+ δ₀*sin.(τ)))
    y = S.R0*(S.ϵ*S.κ*sin.(τ))
    return x,y
end

function (S::MillerShape)(θ)
    δ₀ = asin(S.δ)
    x = S.R0*(1 + S.ϵ * cos(θ + δ₀*sin(θ)))
    y = S.R0*(S.ϵ*S.κ*sin(θ))
    return x,y
end

"""
MillerExtendedHarmonicShape Structure

Defines the Miller Extended Harmonic Plasma Shape Parameterization
> Arbon, Ryan, Jeff Candy, and Emily A. Belli. "Rapidly-convergent flux-surface shape parameterization." Plasma Physics and Controlled Fusion 63.1 (2020): 012001.

Fields:\\
`R0` - Major Radius [m]\\
`Z0` - Elevation [m]
`ϵ`  - Inverse Aspect Ratio a/R0 where a = minor radius\\
`κ`  - Elongation\\
`c0` - Tilt\\
`c`  - Cosine coefficients i.e. acos.([ovality,...])
`s`  - Sine coefficients i.e. asin.([triangularity,-squareness,...])
"""
mutable struct MillerExtendedHarmonicShape{T<:Number} <: PlasmaShape
    R0::T        # Major Radius
    Z0::T        # Elevation
    ϵ::T         # inverse aspect ratio a/R0
    κ::T         # Elongation
    c0::T        # Tilt
    c::Vector{T} # Cosine coefficients acos.([ovality,...])
    s::Vector{T} # Sine coefficients asin.([triangularity,-squareness,...]
end

function MillerExtendedHarmonicShape(R0,ϵ,κ,c::Vector,s::Vector; Z0=zero(R0),c0=zero(R0))
    any(-1.0 .<= s .<= 1.0) && error("Sine coefficients out of range [-1, 1]")
    any(-1.0 .<= c .<= 1.0) && error("Cosine coefficients out of range [-1, 1]")

    MillerExtendedHarmonicShape(promote(R0,Z0,ϵ,κ,c0)...,acos.(c), asin.(s))
end

function MillerExtendedHarmonicShape(R0, ϵ, κ, δ; Z0=zero(R0),tilt=zero(R0), ovality=one(R0), squareness=zero(R0))
    abs(δ) > 1 && error("Triangularity out of range [-1,1]")
    abs(squareness) > 1 && error("Squareness out of range [-1,1]")
    abs(ovality) > 1 && error("Ovality out of range [-1,1]")

    MillerExtendedHarmonicShape(promote(R0,Z0,ϵ,κ,c0)...,acos.([ovality]), asin.([δ,-squareness]))
end

function MillerExtendedHarmonicShape(S::MillerShape)
    s = [S.δ]
    MillerExtendedHarmonicShape(S.R0,zero(S.R0),S.ϵ,S.κ,zero(S.κ),zero(s), asin.(s))
end

elevation(S::MillerExtendedHarmonicShape) = S.Z0
tilt(S::MillerExtendedHarmonicShape) = S.c0
ovality(S::MillerExtendedHarmonicShape) = S.c[1]
triangularity(S::MillerExtendedHarmonicShape) = sin(S.s[1])
squareness(S::MillerExtendedHarmonicShape) = sin(-S.s[2])

function shape(S::MillerExtendedHarmonicShape; N=100)
    x = zeros(N)
    y = zeros(N)
    r = S.ϵ*S.R0
    θ = range(0,2pi,length=N)
    
    for i=1:N
        c_sum = sum(S.c[n]*cos(n*θ[i]) for n=1:length(S.c))
        s_sum = sum(S.s[n]*sin(n*θ[i]) for n=1:length(S.s))
        θ_R = θ[i] + S.c0 + c_sum + s_sum
        x[i] = S.R0 + r*cos(θ_R)
        y[i] = S.Z0 + S.κ*r*sin(θ[i])
    end

    return x,y
end

function (S::MillerExtendedHarmonicShape)(θ)
    r = S.ϵ*S.R0
    c_sum = sum(S.c[n]*cos(n*θ) for n=1:length(S.c))
    s_sum = sum(S.s[n]*sin(n*θ) for n=1:length(S.s))
    θ_R = θ[i] + S.c0 + c_sum + s_sum
    x = S.R0 + r*cos(θ_R)
    y = S.Z0 + S.κ*r*sin(θ)

    return x,y
end

function curvature(S::PlasmaShape,θ)
    xp = ForwardDiff.derivative(x->S(x)[1],θ)
    yp = ForwardDiff.derivative(x->S(x)[2],θ)
    xpp = ForwardDiff.derivative(t->ForwardDiff.derivative(x->S(x)[1],t),θ)
    ypp = ForwardDiff.derivative(t->ForwardDiff.derivative(x->S(x)[2],t),θ)

    κ = abs(yp*xpp - ypp*xp)/(xp^2 + yp^2)^1.5
    return κ
end

