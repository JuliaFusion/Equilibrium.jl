mutable struct EFITEquilibrium{T<:Real, S<:AbstractRange{Float64},
                               R<:AbstractMatrix{Float64},
                               Q<:AbstractVector{Float64}} <: AbstractEquilibrium
    cocos::COCOS       # COCOS
    r::S               # Radius/R range
    z::S               # Elevation/Z range
    psi::S             # Polodial Flux range (polodial flux from magnetic axis)
    psi_rz::R          # Polodial Flux on RZ grid (polodial flux from magnetic axis)
    g::Q               # Polodial Current
    p::Q               # Plasma pressure
    q::Q               # Q profile
    phi::Q             # Electric Potential
    axis::NTuple{2, T} # Magnetic Axis (raxis,zaxis)
    sigma::Int         # sign(dot(J,B))
end

function efit(cc::COCOS, r::S, z::S, psi::S, psi_rz, g, p,
              q, phi, axis::NTuple{2,T}, sigma::Int) where {T,S<:AbstractRange}

    psi_rz_itp = cubic_spline_interpolation((r,z), psi_rz, extrapolation_bc=Flat())
    if step(psi) > 0
        g_itp = cubic_spline_interpolation(psi, g, extrapolation_bc=Flat())
        p_itp = cubic_spline_interpolation(psi, p, extrapolation_bc=Flat())
        q_itp = cubic_spline_interpolation(psi, q, extrapolation_bc=Flat())
        phi_itp = cubic_spline_interpolation(psi, phi, extrapolation_bc=Flat())
    else # cubic_spline_interpolation doesn't like decreasing psi so reverse them
        g_itp = cubic_spline_interpolation(reverse(psi), reverse(g), extrapolation_bc=Flat())
        p_itp = cubic_spline_interpolation(reverse(psi), reverse(p), extrapolation_bc=Flat())
        q_itp = cubic_spline_interpolation(reverse(psi), reverse(q), extrapolation_bc=Flat())
        phi_itp = cubic_spline_interpolation(reverse(psi), reverse(phi), extrapolation_bc=Flat())
    end
    EFITEquilibrium(cc, r, z, psi, psi_rz_itp, g_itp, p_itp, q_itp, phi_itp, axis, Int(sigma))
end

function Base.show(io::IO, N::EFITEquilibrium)
    print(io, "EFITEquilibrium{$(eltype(N.psi))}")
end

Base.broadcastable(N::EFITEquilibrium) = (N,)

cocos(N::EFITEquilibrium) = N.cocos

function (N::EFITEquilibrium)(r,z)
    return N.psi_rz(r,z)
end

function magnetic_axis(N::EFITEquilibrium)
    return N.axis
end

function B0Ip_sign(N::EFITEquilibrium)
    return N.sigma
end

function limits(N::EFITEquilibrium)
    return extrema(N.r), extrema(N.z)
end

function psi_limits(N::EFITEquilibrium)
    return N.psi[1], N.psi[end]
end

function psi_gradient(N::EFITEquilibrium,r,z)
    return SVector{2}(Interpolations.gradient(N.psi_rz, r, z))
end

function poloidal_current(N::EFITEquilibrium, psi)
    return N.g(psi)
end

function poloidal_current_gradient(n::EFITEquilibrium, psi)
    return Interpolations.gradient(n.g, psi)[1]
end

function pressure(N::EFITEquilibrium, psi)
    return N.p(psi)
end

function pressure_gradient(n::EFITEquilibrium, psi)
    return Interpolations.gradient(n.p, psi)[1]
end

function safety_factor(N::EFITEquilibrium, psi)
    return N.q(psi)
end

function electric_potential(N::EFITEquilibrium, psi)
    return N.phi(psi)
end

function electric_potential_gradient(N::EFITEquilibrium, psi)
    return Interpolations.gradient(N.phi, psi)[1]
end
