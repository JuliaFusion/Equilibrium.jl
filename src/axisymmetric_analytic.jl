struct AxisymmetricAnalyticalEquilibrium{T,F<:Function}
    B0::T
    axis::NTuple{2,T}
    sigmaB::Int
    sigmaJ::Int
    q::F
end

function Base.show(io::IO, M::AxisymmetricAnalyticEquilibrium{T,F}) where {T,F}
    print(io, "AxisymmetricAnalyticEquilibrium{$(T)}")
end

Base.broadcastable(M::AxisymmetricAnalyticEquilibrium) = (M,)

function Bfield(M::AxisymmetricAnalyticEquilibrium, r, z)
    Rm = M.axis[1]
    rs = r - Rm
    zs = z - M.axis[2]

    rp = hypot(rs,zs)
    theta = atan(zs,rs)

    A = ((M.B0*Rm)/(Rm - rp*cos(theta)))
    q = M.q(rp)
    B_theta = A*r*M.sigmaJ/(q*Rm)
    s, c = sincos(theta)
    Br = -B_theta*s
    Bz = B_theta*c

    Bt = A*M.sigmaB

    return SVector{3}(Br,Bt,Bz)
end

function Efield(M::AxisymmetricAnalyticEquilibrium, r, z)
    z = zero(r)
    return SVector{3}(z,z,z)
end
