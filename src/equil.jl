mutable struct AxisymmetricEquilibrium{T<:Real, S<:AbstractRange{Float64},
                               R<:AbstractArray{Float64,2},
                               Q<:AbstractArray{Float64,1}}
    r::S               # Radius/R range
    z::S               # Elevation/Z range
    psi::S             # Polodial Flux range (polodial flux from magnetic axis)
    psi_rz::R          # Polodial Flux on RZ grid (polodial flux from magnetic axis)
    g::Q               # Polodial Current
    p::Q               # Plasma pressure
    q::Q               # Q profile
    phi::Q             # Electric Potential
    b::R               # Magnetic field magnitude
    j::R               # Plasma Current magnitude
    axis::NTuple{2, T} # Magnetic Axis (raxis,zaxis)
    sigma::Int         # sign(dot(J,B))
    alpha::Int         # Poloidal Flux type "Disk: 1" or "Ribbon: -1"
    flux::T            # Enclosed Poloidal Flux
end

function AxisymmetricEquilibrium(r::AbstractRange{T}, z::AbstractRange{T}, psi::AbstractRange{T}, psi_rz, g, p, q, phi, axis::NTuple{2,T}, psi_mag, psi_bdry) where {T <: Real}

    flux = psi_bdry - psi_mag
    if iszero(psi_mag)
        alpha = -1
    else
        alpha = 1
    end

    psi_rz_itp = CubicSplineInterpolation((r,z), psi_rz, extrapolation_bc=Flat())
    g_itp = CubicSplineInterpolation(psi, g, extrapolation_bc=Flat())
    p_itp = CubicSplineInterpolation(psi, p, extrapolation_bc=Flat())
    q_itp = CubicSplineInterpolation(psi, q, extrapolation_bc=Flat())
    phi_itp = CubicSplineInterpolation(psi, phi, extrapolation_bc=Flat())

    b = [norm(Bfield(psi_rz_itp,g_itp,rr,zz,alpha)) for rr in r, zz in z]
    b_itp = CubicSplineInterpolation((r,z),b,extrapolation_bc=Flat())

    j = [norm(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz,alpha)) for rr in r, zz in z]
    j_itp = CubicSplineInterpolation((r,z),j,extrapolation_bc=Flat())

    rr = axis[1] + (r[end] - axis[1])/10
    zz = axis[2] + (z[end] - axis[2])/10

    sigma = Int(sign(dot(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz,alpha), Bfield(psi_rz_itp,g_itp,rr,zz,alpha))))

    AxisymmetricEquilibrium(r, z, psi, psi_rz_itp, g_itp, p_itp, q_itp, phi_itp, b_itp, j_itp, axis, sigma, alpha, flux)
end

function Base.show(io::IO, M::AxisymmetricEquilibrium)
    print(io, "AxisymmetricEquilibrium{$(eltype(M.psi))}")
end

Base.broadcastable(M::AxisymmetricEquilibrium) = (M,)

struct EMFields{S<:Number,T<:AbstractVector}
    psi::S
    g::S
    B::T
    E::T
end

function Base.show(io::IO, EM::EMFields)
    print(io, "EMFields")
    print(io, " B = $(round.(EM.B,digits=3)) [T]")
    print(io, " E = $(round.(EM.E,digits=3)) [V/m]")
end

function EMFields(psi_rz, g, phi, r, z, alpha)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = alpha*SVector{2}(Interpolations.gradient(psi_rz, r, z))
    grad_psi_norm = norm(grad_psi)

    # Calculate B-Field
    BR = -grad_psi[2]/r
    Bz = grad_psi[1]/r
    Bt = gval/r
    Bpol = sqrt(BR^2 + Bz^2)

    # Calculate E-Field
    Er = alpha*r*Bpol*Interpolations.gradient(phi, psi)[1]
    ER = -Er*grad_psi[1]/grad_psi_norm # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = -Er*grad_psi[2]/grad_psi_norm # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = zero(Ez)

    return EMFields(psi, gval, SVector{3}(BR,Bt,Bz), SVector{3}(ER,Et,Ez))
end

function EMFields(M::AxisymmetricEquilibrium, r, z)
    EMFields(M.psi_rz, M.g, M.phi, r, z, M.alpha)
end

function fields(M::AxisymmetricEquilibrium, r, z)
    F = EMFields(M.psi_rz, M.g, M.phi, r, z, M.alpha)
    return F
end

function fields(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    F = fields(M,r,z)
    sp, cp = sincos(phi)
    B = SVector{3}(F.B[1]*cp - F.B[2]*sp,F.B[1]*sp + F.B[2]*cp,F.B[3])
    E = SVector{3}(F.E[1]*cp - F.E[2]*sp,F.E[1]*sp + F.E[2]*cp,F.E[3])
    return EMFields(F.psi, F.g, B, E)
end

function Bfield(psi_rz, g, r, z, alpha)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = alpha*SVector{2}(Interpolations.gradient(psi_rz, r, z))

    br = -grad_psi[2]/r
    bz = grad_psi[1]/r
    bt = gval/r

    return SVector{3}(br,bt,bz)
end

function Bfield(M::AxisymmetricEquilibrium, r, z)
    B = Bfield(M.psi_rz, M.g, r, z, M.alpha)
    return B
end

function Bfield(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    B = Bfield(M,r,z)
    sp, cp = sincos(phi)
    B_xyz = SVector{3}(B[1]*cp - B[2]*sp, B[1]*sp + B[2]*cp, B[3])
    return B_xyz
end

function Jfield(psi_rz, g, p, r, z, alpha)
    psi = psi_rz(r,z)
    gval = g(psi)
    grad_psi = alpha*SVector{2}(Interpolations.gradient(psi_rz, r, z))

    gp = alpha*Interpolations.gradient(g, psi)[1]
    pp = alpha*Interpolations.gradient(p, psi)[1]

    jr = -gp*grad_psi[2]/(r*mu0)
    jz = gp*grad_psi[1]/(r*mu0)
    jt = r*pp + gval*gp/(r*mu0)

    return SVector{3}(jr,jt,jz)
end

function Jfield(M::AxisymmetricEquilibrium, r, z)
    J = Jfield(M.psi_rz, M.g, M.p, r, z, M.alpha)
    return J
end

function Jfield(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    J = Jfield(M,r,z)
    sp, cp = sincos(phi)
    J_xyz = SVector{3}(J[1]*cp - J[2]*sp, J[1]*sp + J[2]*cp, J[3])
    return J_xyz
end

function Efield(psi_rz, phi, r, z, Bpol, alpha)
    psi = psi_rz(r,z)
    grad_psi = alpha*SVector{2}(Interpolations.gradient(psi_rz, r, z))
    grad_psi = grad_psi/norm(grad_psi)

    Er = alpha*r*Bpol*Interpolations.gradient(phi, psi)[1]

    ER = -Er*grad_psi[1] # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = -Er*grad_psi[2] # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = zero(Ez)

    return SVector{3}(ER, Et, Ez)
end

function Efield(M::AxisymmetricEquilibrium, r, z)
    B = Bfield(M, r, z)
    Bpol = sqrt(B[1]^2 + B[3]^2)
    E = Efield(M.psi_rz, M.phi, r, z, Bpol, M.alpha)
    return E
end

function Efield(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    E = Efield(M,r,z)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1]*cp - E[2]*sp, E[1]*sp + E[2]*cp, E[3])
    return E_xyz
end

function Efield(M::AxisymmetricEquilibrium, r, z, vrot::AbstractVector)
    psi = M.psi_rz(r,z)
    gval = M.g(psi)
    grad_psi = M.alpha*SVector{2}(Interpolations.gradient(M.psi_rz, r, z))
    B = SVector{3}(-grad_psi[2]/r, grad_psi[1]/r, gval/r)

    grad_psi = grad_psi/norm(grad_psi)
    qval = M.q(psi)
    dp = M.alpha*Interpolations.gradient(M.p, psi)[1]
    gradp = SVector{3}(dp*grad_psi[1], zero(grad_psi[1]), dp*grad_psi[2])

    E = cross(vrot,B) .+ gradp/qval
    return E
end

function Efield(M::AxisymmetricEquilibrium, x, y, z, vrot::AbstractVector)
    r = hypot(x,y)
    phi = atan(y,x)
    E = Efield(M,r,z,vrot)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1]*cp - E[2]*sp, E[1]*sp + E[2]*cp, E[3])
    return E_xyz
end

function rho_p(M::AxisymmetricEquilibrium, r, z)
    sqrt((M.psi_rz(r,z) - M.psi[1])/M.flux)
end

function rho_p(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    rho_p(M,r,z)
end

function gradB(M::AxisymmetricEquilibrium, r, z)
    gB_rz = Interpolations.gradient(M.b,r,z)
    SVector{3}(gB_rz[1],0.0,gB_rz[2])
end

function gradB(M::AxisymmetricEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    gB_rz = Interpolations.gradient(M.b,r,z)
    gB = SVector{3}(gB_rz[1],0.0,gB_rz[2])
    sp, cp = sincos(phi)
    gB_xyz = SVector{3}(gB[1]*cp - gB[2]*sp, gB[1]*sp + gB[2]*cp, gB[3])
end

function curlB(M::AxisymmetricEquilibrium, r, z)
    J = ForwardDiff.jacobian(x->Bfield(M,x[1],x[3]),SVector{3}(r,0.0,z))
    B = Bfield(M,r,z)
    return SVector{3}(J[3,2]/r - J[2,3], J[1,3] - J[3,1], (B[2]/r + J[2,1]) - J[1,2])
end

function curlB(M::AxisymmetricEquilibrium, x, y, z)
    J = ForwardDiff.jacobian(x->Bfield(M,x[1],x[2],x[3]),SVector{3}(x,y,z))
    return SVector{3}(J[3,2] - J[2,3], J[1,3]-J[3,1], J[2,1] - J[1,2])
end

