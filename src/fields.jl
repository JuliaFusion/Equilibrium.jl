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

function EMFields(M::T, r, z) where T<:AbstractEquilibrium
    psi = M(r,z)
    gval = poloidal_current(M,r,z)
    B = Bfield(M, r, z)
    E = Efield(M, r, z)
    return EMFields(psi, gval, B, E)
end

function fields(M::AbstractEquilibrium, r, z)
    return EMFields(M, r, z)
end

function fields(M::AbstractEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    F = fields(M,r,z)
    sp, cp = sincos(phi)
    B = SVector{3}(F.B[1]*cp - F.B[2]*sp,F.B[1]*sp + F.B[2]*cp,F.B[3])
    E = SVector{3}(F.E[1]*cp - F.E[2]*sp,F.E[1]*sp + F.E[2]*cp,F.E[3])
    return EMFields(F.psi, F.g, B, E)
end

function Bfield(M::T, r, z) where T<:AbstractEquilibrium
    psi = M(r,z)
    gval = poloidal_current(M,psi)
    grad_psi = psi_gradient(M,r,z)
    grad_psi_norm = norm(grad_psi)

    # Calculate B-Field
    cc = cocos(M)
    cocos_factor = cc.sigma_RpZ*cc.sigma_Bp/((2pi)^cc.exp_Bp)
    BR =  cocos_factor*grad_psi[2]/r
    Bz = -cocos_factor*grad_psi[1]/r
    Bt = gval/r

    B = SVector{3}(cylindrical_cocos(cc, BR, Bt, Bz))
    return B
end

function Bfield(M::T, x, y, z) where T<:AbstractEquilibrium
    r = hypot(x,y)
    phi = atan(y,x)
    B = Bfield(M,r,z)
    sp, cp = sincos(phi)
    B_xyz = SVector{3}(B[1]*cp - B[2]*sp, B[1]*sp + B[2]*cp, B[3])
    return B_xyz
end

function poloidal_Bfield(M::T, r, z) where T<:AbstractEquilibrium
    rmaxis, zmaxis = magnetic_axis(M)
    cc = cocos(M)
    B = Bfield(M, r, z)

    ir, iphi, iz = cylindrical_cocos_indices(cc)
    sign_theta = cc.sigma_RpZ*cc.sigma_rhotp # + CW, - CCW
    sign_Bp    = sign_theta*sign((z-zmaxis)*B[ir] - (r-rmaxis)*B[iz]) # sign(theta)*sign(r x B)
    Bpol = sign_Bp*sqrt(B[ir]^2 + B[iz]^2)
    return Bpol
end

function Jfield(M::T, r, z; curl=false) where T<:AbstractEquilibrium
    if curl
        return curlB(M, r, z)/mu0
    end

    psi = M(r,z)
    gval = poloidal_current(M,psi)
    grad_psi = psi_gradient(M,r,z)

    gp = poloidal_current_gradient(M,psi)
    pp = pressure_gradient(M,psi)

    cc = cocos(M)

    cocos_factor = cc.sigma_RpZ
    jr = -cocos_factor*gp*grad_psi[2]/(r*mu0)
    jz =  cocos_factor*gp*grad_psi[1]/(r*mu0)

    cocos_factor = -cc.sigma_Bp*((2pi)^cc.exp_Bp)
    jt = cocos_factor*(r*pp + gval*gp/(r*mu0))

    return SVector{3}(cylindrical_cocos(cc,jr,jt,jz))
end

function Jfield(M::T, x, y, z; kwargs...) where T<:AbstractEquilibrium
    r = hypot(x,y)
    phi = atan(y,x)
    J = Jfield(M,r,z; kwargs...)
    sp, cp = sincos(phi)
    J_xyz = SVector{3}(J[1]*cp - J[2]*sp, J[1]*sp + J[2]*cp, J[3])
    return J_xyz
end

function poloidal_Jfield(M::T, r, z) where T<:AbstractEquilibrium
    rmaxis, zmaxis = magnetic_axis(M)
    cc = cocos(M)
    B = Jfield(M, r, z)

    ir, iphi, iz = cylindrical_cocos_indices(cc)
    sign_theta = cc.sigma_RpZ*cc.sigma_rhotp # + CW, - CCW
    sign_Jp    = sign_theta*sign((z-zmaxis)*J[ir] - (r-rmaxis)*J[iz]) # sign(theta)*sign(r x B)
    Jpol = sign_Jp*sqrt(J[ir]^2 + J[iz]^2)
    return Jpol
end

function Efield(M::T, r, z) where T<:AbstractEquilibrium

    psi = M(r,z)
    phi_grad = phi_gradient(M, psi)
    if phi_grad == 0
        return SVector{3}(zero(r),zero(r),zero(r))
    end

    Bpol = poloidal_Bfield(M, r, z) #signed

    grad_psi = psi_gradient(M,r,z)
    grad_psi = grad_psi/norm(grad_psi)

    # is Bpol signed in this equation?
    Er = -r*Bpol*phi_grad

    ER = Er*grad_psi[1] # Er*dpsi/dR = (-dphi/dpsi)*(dpsi/dR)
    Ez = Er*grad_psi[2] # Er*dpsi/dz = (-dphi/dpsi)*(dpsi/dz)
    Et = zero(Ez)

    return SVector{3}(cylindrical_cocos(cc, ER, Et, Ez))
end

function Efield(M::T, x, y, z) where T<:AbstractEquilibrium
    r = hypot(x,y)
    phi = atan(y,x)
    E = Efield(M,r,z)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1]*cp - E[2]*sp, E[1]*sp + E[2]*cp, E[3])
    return E_xyz
end

function Efield(M::T, r, z, vrot::AbstractVector) where T<:AbstractEquilibrium
    cc = cocos(M)
    psi = M(r,z)
    gval = poloidal_current(M,psi)
    B = Bfield(M, r, z)

    grad_psi = psi_gradient(M,r,z)
    grad_psi = grad_psi/norm(grad_psi)

    qval = safety_factor(M,psi)
    dp = pressure_gradient(M, psi)

    gradp = SVector{3}(cylindrical_cocos(cc, dp*grad_psi[1], zero(grad_psi[1]), dp*grad_psi[2]))

    E = cross(vrot,B) .+ gradp/qval

    return E
end

function Efield(M::T, x, y, z, vrot::AbstractVector) where T<:AbstractEquilibrium
    r = hypot(x,y)
    phi = atan(y,x)
    E = Efield(M,r,z,vrot)
    sp, cp = sincos(phi)
    E_xyz = SVector{3}(E[1]*cp - E[2]*sp, E[1]*sp + E[2]*cp, E[3])
    return E_xyz
end

function rho_p(M::T, psi) where T<:AbstractEquilibrium
    psimag,psibry = psi_limits(M)
    flux = abs(psibry - psimag)
    return sqrt((psi - psimag)/flux)
end

function rho_p(M::T, r, z) where T<:AbstractEquilibrium
    rho_p(M,M(r,z))
end

function rho_p(M::T, x, y, z) where T<:AbstractEquilibrium
    r = hypot(x,y)
    rho_p(M,r,z)
end

function toroidal_flux(M::T,psi) where T<:AbstractEquilibrium
    cc = cocos(M)
    cocos_factor = cc.sigma_Bp*cc.sigma_rhotp*(2pi)^(1 - cc.exp_Bp)

    psimag, psibry = psi_limits(M)
    phi = cocos_factor*hquadrature(x->safety_factor(M,x), psimag, psi)[1]
    return phi
end

function safety_factor(M::AbstractEquilibrium, psi)
    g = poloidal_current(M,psi)
    fs = boundary(M, psi)
    q = (g/(2pi))*average(fs,(x,y)->inv(poloidal_Bfield(M,x,y)*x^2))*circumference(fs)
    return q
end

function gradB(M::AbstractEquilibrium, r, z)
    gB_rz = ForwardDiff.gradient(x->norm(Bfield(M,x[1],x[2])), SVector{2}(r,z))
    return SVector{3}(cylindrical_cocos(cocos(M), gB_rz[1], 0.0, gB_rz[2]))
end

function gradB(M::AbstractEquilibrium, x, y, z)
    r = hypot(x,y)
    phi = atan(y,x)
    gB = gradB(M, r, z)
    sp, cp = sincos(phi)
    gB_xyz = SVector{3}(gB[1]*cp - gB[2]*sp, gB[1]*sp + gB[2]*cp, gB[3])
end

function curlB(M::AbstractEquilibrium, r, z)
    cc = cocos(M)
    J = ForwardDiff.jacobian(x->Bfield(M,x[1],x[3]),SVector{3}(r,0.0,z))
    B = Bfield(M,r,z)
    return cc.sigma_RpZ*SVector{3}(cylindrical_cocos(cc, J[3,2]/r - J[2,3], J[1,3] - J[3,1], (B[2]/r + J[2,1]) - J[1,2]/r))
end

function curlB(M::AbstractEquilibrium, x, y, z)
    J = ForwardDiff.jacobian(x->Bfield(M,x[1],x[2],x[3]),SVector{3}(x,y,z))
    return SVector{3}(J[3,2] - J[2,3], J[1,3]-J[3,1], J[2,1] - J[1,2])
end
