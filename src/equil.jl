mutable struct AxisymmetricEquilibrium{T<:Real, S<:AbstractRange{Float64},
                               R<:AbstractArray{Float64,2},
                               Q<:AbstractArray{Float64,1}}
    r::S               # Radius/R range
    z::S               # Elevation/Z range
    psi::S             # "Ribbon" Polodial Flux range (polodial flux from magnetic axis)
    psi_rz::R          # "Ribbon" Polodial Flux on RZ grid (polodial flux from magnetic axis)
    g::Q               # Polodial Current
    p::Q               # Plasma pressure
    q::Q               # Q profile
    phi::Q             # Electric Potential
    b::R               # Magnetic field magnitude
    j::R               # Plasma Current magnitude
    axis::NTuple{2, T} # Magnetic Axis (raxis,zaxis)
    sigma::Int         # sign(dot(J,B))
    flux::T            # Enclosed Poloidal Flux
end

function AxisymmetricEquilibrium(r::AbstractRange{T}, z::AbstractRange{T}, psi::AbstractRange{T}, psi_rz, g, p, q, phi, axis::NTuple{2,T}, flux) where {T <: Real}

    psi_rz_itp = CubicSplineInterpolation((r,z), psi_rz, extrapolation_bc=Flat())
    g_itp = CubicSplineInterpolation(psi, g, extrapolation_bc=Flat())
    p_itp = CubicSplineInterpolation(psi, p, extrapolation_bc=Flat())
    q_itp = CubicSplineInterpolation(psi, q, extrapolation_bc=Flat())
    phi_itp = CubicSplineInterpolation(psi, phi, extrapolation_bc=Flat())

    b = [norm(Bfield(psi_rz_itp,g_itp,rr,zz)) for rr in r, zz in z]
    b_itp = CubicSplineInterpolation((r,z),b,extrapolation_bc=Flat())

    j = [norm(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz)) for rr in r, zz in z]
    j_itp = CubicSplineInterpolation((r,z),j,extrapolation_bc=Flat())

    rr = axis[1] + (r[end] - axis[1])/10
    zz = axis[2] + (z[end] - axis[2])/10

    sigma = Int(sign(dot(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz), Bfield(psi_rz_itp,g_itp,rr,zz))))

    AxisymmetricEquilibrium(r, z, psi, psi_rz_itp, g_itp, p_itp, q_itp, phi_itp, b_itp, j_itp, axis, sigma, flux)
end

function Base.show(io::IO, M::AxisymmetricEquilibrium)
    print(io, "AxisymmetricEquilibrium{$(eltype(M.psi))}")
end

struct EMFields{T<:AbstractVector}
    B::T
    E::T
end

function Base.show(io::IO, EM::EMFields)
    print(io, "EMFields")
    print(io, " B = $(round.(EM.B,digits=3)) [T]")
    print(io, " E = $(round.(EM.E,digits=3)) [V/m]")
end

function EMFields(psi_rz, g, phi, r, z)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = Interpolations.gradient(psi_rz, r, z)
    grad_psi_norm = norm(grad_psi)

    # Calculate B-Field
    BR = grad_psi[2]/r
    Bz = -grad_psi[1]/r
    Bt = gval/r
    Bpol = sqrt(BR^2 + Bz^2)

    # Calculate E-Field
    Er = -r*Bpol*Interpolations.gradient(phi, psi)[1]
    ER = Er*grad_psi[1]/grad_psi_norm # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = Er*grad_psi[2]/grad_psi_norm # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = 0.0

    return EMFields([BR,Bt,Bz],[ER,Et,Ez])
end

function EMFields(M::AxisymmetricEquilibrium, r, z)
    EMFields(M.psi_rz, M.g, M.phi, r, z)
end

function fields(M::AxisymmetricEquilibrium, r, z)
    EMFields(M.psi_rz, M.g, M.phi, r, z)
end

function Bfield(psi_rz, g, r, z)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = Interpolations.gradient(psi_rz, r, z)

    br = grad_psi[2]/r
    bz = -grad_psi[1]/r
    bt = gval/r

    return [br,bt,bz]
end

function Bfield(M::AxisymmetricEquilibrium, r, z)
    Bfield(M.psi_rz, M.g, r, z)
end

function Jfield(psi_rz, g, p, r, z)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = Interpolations.gradient(psi_rz, r, z)

    gp = -Interpolations.gradient(g, psi)[1]
    pp = -Interpolations.gradient(p, psi)[1]

    jr = gp*grad_psi[2]/(r*mu0)
    jz = -gp*grad_psi[1]/(r*mu0)
    jt = r*pp + gval*gp/(r*mu0)

    return [jr,jt,jz]
end

function Jfield(M::AxisymmetricEquilibrium, r, z)
    Jfield(M.psi_rz, M.g, M.p, r, z)
end

function Efield(psi_rz, phi, r, z, Bpol)
    psi = psi_rz(r,z)
    grad_psi = Interpolations.gradient(psi_rz, r, z)
    grad_psi = grad_psi/norm(grad_psi)

    Er = -r*Bpol*Interpolations.gradient(phi, psi)[1]

    ER = Er*grad_psi[1] # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = Er*grad_psi[2] # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = 0.0

    return [ER, Et, Ez]
end

function Efield(M::AxisymmetricEquilibrium, r, z)
    B = Bfield(M, r, z)
    Bpol = sqrt(B[1]^2 + B[3]^2)
    Efield(M.psi_rz, M.phi, r, z, Bpol)
end
