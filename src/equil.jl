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

    psi_max = maximum(psi_rz)
    dpsi = step(psi)
    psi_ext = psi[1]:dpsi:psi_max
    psi_ext = range(psi_ext[1],stop=psi_ext[end],length=length(psi_ext))
    g_ext = [i <= length(g) ? g[i] : g[end] for i=1:length(psi_ext)]
    p_ext = [i <= length(p) ? p[i] : 0.0 for i=1:length(psi_ext)]
    q_ext = [i <= length(q) ? q[i] : q[end] for i=1:length(psi_ext)]
    phi_ext = [i <= length(phi) ? phi[i] : phi[end] for i=1:length(psi_ext)]

    psi_rz_itp = scale(interpolate(psi_rz, BSpline(Cubic(Line())), OnGrid()), r, z)
    g_itp = scale(interpolate(g_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)
    p_itp = scale(interpolate(p_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)
    q_itp = scale(interpolate(q_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)
    phi_itp = scale(interpolate(phi_ext, BSpline(Cubic(Line())), OnGrid()), psi_ext)

    b = [norm(Bfield(psi_rz_itp,g_itp,rr,zz)) for rr in r, zz in z]
    b_itp = scale(interpolate(b, BSpline(Cubic(Line())), OnGrid()), r, z)

    j = [norm(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz)) for rr in r, zz in z]
    j_itp = scale(interpolate(j, BSpline(Cubic(Line())), OnGrid()), r, z)

    rr = axis[1] + (r[end] - axis[1])/10
    zz = axis[2] + (z[end] - axis[2])/10

    sigma = Int(sign(dot(Jfield(psi_rz_itp,g_itp,p_itp,rr,zz), Bfield(psi_rz_itp,g_itp,rr,zz))))

    AxisymmetricEquilibrium(r, z, psi_ext, psi_rz_itp, g_itp, p_itp, q_itp, phi_itp, b_itp, j_itp, axis, sigma, flux)
end

function Base.show(io::IO, M::AxisymmetricEquilibrium)
    print(io, "AxisymmetricEquilibrium{$(eltype(M.psi))}")
end

function Bfield(psi_rz, g, r, z)
    psi = psi_rz[r,z]
    gval = g[psi]
    grad_psi = gradient(psi_rz, r, z)

    br = grad_psi[2]/r
    bz = -grad_psi[1]/r
    bt = gval/r

    return [br,bt,bz]
end

function Bfield(M::AxisymmetricEquilibrium, r, z)
    Bfield(M.psi_rz, M.g, r, z)
end

function Jfield(psi_rz, g, p, r, z)
    psi = psi_rz[r,z]
    gval = g[psi]
    grad_psi = gradient(psi_rz, r, z)

    gp = -gradient(g, psi)[1]
    pp = -gradient(p, psi)[1]

    jr = gp*grad_psi[2]/(r*mu0)
    jz = -gp*grad_psi[1]/(r*mu0)
    jt = r*pp + gval*gp/(r*mu0)

    return [jr,jt,jz]
end

function Jfield(M::AxisymmetricEquilibrium, r, z)
    Jfield(M.psi_rz, M.g, M.p, r, z)
end

function Efield(psi_rz, phi, r, z, Bpol)
    psi = psi_rz[r,z]
    grad_psi = gradient(psi_rz, r, z)
    grad_psi = grad_psi/norm(grad_psi)

    Er = -r*Bpol*gradient(phi, psi)[1]

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
