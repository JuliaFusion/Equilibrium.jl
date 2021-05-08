"""
COCOS Structure

`cocos::Int`           = COCOS ID number\\
`exp_Bp::Int`          = 0 or 1, depending if psi is already divided by 2pi or not, respectively\\
`sigma_Bp::Int`        = +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive\\
`sigma_RpZ::Int`       = +1 or -1, depending if (R,phi,Z) is right-handed or (R,Z,phi), respectively\\
`sigma_rhotp::Int`     = +1 or -1, depending if (rho, theta, phi) is right-handed or (rho,phi,theta), repectively\\
`sign_q_pos::Int`      = +1 or -1, depending if q is positive or negative with Ip and B0 positive\\
`sign_pprime_pos::Int` = +1 or -1, depending if dp/dpsi is positive or negative with Ip and B0 positive
"""
struct COCOS
    cocos::Int           # COCOS ID number
    exp_Bp::Int          # 0 or 1, depending if psi is already divided by 2pi or not, respectively
    sigma_Bp::Int        # +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive
    sigma_RpZ::Int       # +1 or -1, depending if (R,phi,Z) is right-handed or (R,Z,phi), respectively
    sigma_rhotp::Int     # +1 or -1, depending if (rho, theta, phi) is right-handed or (rho,phi,theta), repectively
    sign_q_pos::Int      # +1 or -1, depending if q is positive or negative with Ip and B0 positive
    sign_pprime_pos::Int # +1 or -1, depending if dp/dpsi is positive or negative with Ip and B0 positive
end

cocos(CC::COCOS) = CC

function Base.show(io::IO, CC::COCOS)
    println(io, "COCOS = $(CC.cocos)")
    println(io, " e_Bp  = $(CC.exp_Bp)")
    println(io, " σ_Bp  = $(CC.sigma_Bp)")

    rpz = Dict(1=>"(R,Φ,Z)", -1=>"(R,Z,Φ)")
    rpz_dir = Dict(1=>"CCW", -1=>"CW")
    rhotp = Dict(1=>"(ρ,θ,Φ)", -1=>"(ρ,Φ,θ)")
    rhotp_dir = Dict(1=>"CW", -1=>"CCW")

    println(io, " σ_RΦZ = $(rpz[CC.sigma_RpZ]): $(CC.sigma_RpZ)")
    println(io, " σ_ρθΦ = $(rhotp[CC.sigma_rhotp]): $(CC.sigma_rhotp)")
    println(io, " Φ from top: $(rpz_dir[CC.sigma_RpZ])")
    println(io, " θ from front: $(rhotp_dir[CC.sigma_RpZ*CC.sigma_rhotp])")

    inc = Dict(1=>"Increasing", -1=>"Decreasing")
    println(io, " ψ_ref: $(inc[CC.sigma_Bp]) assuming +Ip, +B0")
    println(io, " sign(q) = $(CC.sign_q_pos) assuming +Ip, +B0")
    print(io, " sign(p') = $(CC.sign_pprime_pos) assuming +Ip, +B0")
end

Base.broadcastable(CC::COCOS) = (CC,)

function _cylindrical_cocos(::Val{1}, r, phi, z)
    return (r, phi, z)
end

function _cylindrical_cocos(::Val{-1}, r, phi, z)
    return (r, z, phi)
end

"""
    cylindrical_cocos(c::COCOS, r, phi, z) -> NTuple{3}

Returns the right-handed cylindrical coordinate according to the provided COCOS.

`cocos.sigma_RpZ = +1 -> (r, phi, z)`\\
`cocos.sigma_RpZ = -1 -> (r, z, phi)`

"""
function cylindrical_cocos(c::COCOS, r, phi, z)
    return _cylindrical_cocos(Val(c.sigma_RpZ), r, phi, z)
end

function _cylindrical_cocos_indices(::Val{1})
    return (1,2,3)
end

function _cylindrical_cocos_indices(::Val{-1})
    return (1,3,2)
end

"""
    cylindrical_cocos_indices(c::COCOS) -> NTuple{3}

Returns the indices of the r, phi, and z coordinates, respectively, relative to the provided COCOS.

`cocos.sigma_RpZ = +1 -> (1,2,3)`\\
`cocos.sigma_RpZ = -1 -> (1,3,2)`

"""
function cylindrical_cocos_indices(c::COCOS)
    return _cylindrical_cocos_indices(Val(c.sigma_RpZ))
end

function _poloidal_cocos(::Val{1}, rho, theta, phi)
    return (rho, theta, phi)
end

function _poloidal_cocos(::Val{-1}, rho, theta, phi)
    return (rho, phi, theta)
end

"""
    poloidal_cocos(c::COCOS, rho, theta, phi) -> NTuple{3}

Returns the right-handed poloidal coordinate according to the provided COCOS.

cocos.sigma_rhotp = +1 -> (rho, theta, phi)\\
cocos.sigma_rhotp = -1 -> (rho, phi, theta)
"""
function poloidal_cocos(c::COCOS, rho, theta, phi)
    return _poloidal_cocos(Val(c.sigma_rhotp), rho, theta, phi)
end

function _poloidal_cocos_indices(::Val{1})
    return (1,2,3)
end

function _poloidal_cocos_indices(::Val{-1})
    return (1, 3, 2)
end

"""
    poloidal_cocos_indices(c::COCOS) -> NTuple{3}

Returns the indices of the rho, theta, and phi coordinates, respectively, relative to the provided COCOS.

cocos.sigma_rhotp = +1 -> (1,2,3)\\
cocos.sigma_rhotp = -1 -> (1,3,2)
"""
function poloidal_cocos_indices(c::COCOS, rho, theta, phi)
    return _poloidal_cocos_indices(Val(c.sigma_rhotp), rho, theta, phi)
end

"""
    cocos(cocos_ID) -> COCOS

Returns `COCOS` structure given the `cocos_ID` number
"""
function cocos(cocos_in)

    exp_Bp = cocos_in >= 11 ? 1 : 0

    if cocos_in ∈ (1,11)
        # ITER, Boozer are COCOS=11
        return COCOS(cocos_in,exp_Bp,1,1,1,1,-1)
    elseif cocos_in ∈ (2,12)
        # CHEASE, ONETWO, Hinton-Hazeltine, LION is COCOS=2
        return COCOS(cocos_in,exp_Bp,1,-1,1,1,-1)
    elseif cocos_in ∈ (3,13)
        # Freidberg, CAXE, KINX, EFIT are COCOS=3
        # EU-ITM up to end of 2011 is COCOS=13
        return COCOS(cocos_in,exp_Bp,-1,1,-1,-1,1)
    elseif cocos_in ∈ (4,14)
        return COCOS(cocos_in,exp_Bp,-1,-1,-1,-1,1)
    elseif cocos_in ∈ (5,15)
        return COCOS(cocos_in,exp_Bp,1,1,-1,-1,-1)
    elseif cocos_in ∈ (6,16)
        return COCOS(cocos_in,exp_Bp,1,-1,-1,-1,-1)
    elseif cocos_in ∈ (7,17)
        return COCOS(cocos_in,exp_Bp,-1,1,1,1,1)
    elseif cocos_in ∈ (8,18)
        return COCOS(cocos_in,exp_Bp,-1,-1,1,1,1)
    else
        throw(ArgumentError("COCOS = $cocos_in does not exist"))
    end
end

"""
    check_cocos(sigma_B0, sigma_Ip, sigma_F, sigma_pprime, sigma_q, sigma_dpsi, cc::COCOS; verbose = false) -> Bool

Returns True if equilibrium quantities are consistant with given COCOS

Arguments:\\
`sigma_B0` - Toroidal magnetic field sign\\
`sigma_Ip` - Plasma current sign\\
`sigma_F` - Poloidal current sign\\
`sigma_pprime` - Pressure gradient sign\\
`sigma_dpsi` - dpsi/dr sign\\
`cc::Union{Int,COCOS}` - COCOS structure or ID
"""
function check_cocos(sigma_B0, sigma_Ip, sigma_F, sigma_pprime,
                     sigma_q, sigma_dpsi,
                     cc::Union{Int,COCOS}; verbose=false)

    cc = cocos(cc)

    valid = true
    if sigma_q*cc.sigma_rhotp*sigma_Ip*sigma_B0 < 0
        verbose && @warn "sigma_q($sigma_q) ≠ sigma_rhotp($(cc.sigma_rhotp))*sigma_Ip($sigma_Ip)*sigma_B0($sigma_B0)"
        valid = false
    end

    if sigma_F*sigma_B0 < 0
        verbose && @warn "Signs of F and B0 are not consistant"
        valid = false
    end

    if sigma_dpsi*cc.sigma_Bp*sigma_Ip < 0
        if sigma_dpsi > 0
            verbose && @warn "psi should be decreasing with sign(Ip) = $(sigma_Ip) for COCOS = $(cc.cocos)"
        else
            verbose && @warn "psi should be increasing with sign(Ip) = $(sigma_Ip) for COCOS = $(cc.cocos)"
        end
        valid = false
    elseif sigma_pprime*sigma_Ip*cc.sigma_Bp > 0
        verbose && @warn "sign(pprime) should be $(-sigma_Ip*cc.sigma_Bp)"
        valid = false
    end

    return valid
end

"""
    check_cocos(B0, Ip, F::AbstractVector, pprime::AbstractVector, q::AbstractVector, psi::AbstracVector, cc::COCOS; verbose = false) -> Bool

Returns True if equilibrium quantities are consistant with given COCOS

Arguments:\\
`B0` - Toroidal magnetic field\\
`Ip` - Plasma current\\
`F::AbstractVector` - Poloidal current as a function of `psi`\\
`pprime::AbstracVector` - Pressure gradient w.r.t. `psi` as a function of `psi`\\
`psi::AbstractVector` - Poloidal flux\\
`cc::Union{Int,COCOS}` - COCOS structure or ID
"""
function check_cocos(B0, Ip, F::AbstractVector, pprime::AbstractVector,
                     q::AbstractVector, psi::AbstractVector,
                     cc::Union{Int,COCOS}; kwargs...)

    sigma_B0 = sign(B0)
    sigma_Ip = sign(Ip)
    sigma_F = sign(sum(F))
    sigma_pprime = sign(sum(pprime))
    sigma_q = sign(q[end])
    sigma_dpsi = sign(psi[end] - psi[1])

    return check_cocos(sigma_B0, sigma_Ip, sigma_F, sigma_pprime, sigma_q, sigma_dpsi, cc; kwargs...)

end

"""
    transform_cocos(cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS};
                    sigma_Ip = nothing, sigma_B0=nothing,
                    ld = (1,1), lB = (1,1), exp_mu0 = (0,0)) -> Dict

Returns a dictionary of the multiplicative factors to transform COCOS from `cc_in` to `cc_out`

Arguments:\\
`sigma_Ip::Union{NTuple{2,Int},Nothing}` - A tuple of the (Input, Output) current sign or nothing\\
`sigma_B0::Union{NTuple{2,Int},Nothing}` - A tuple of the (Input, Output) toroidal field sign or nothing\\
`ld::NTuple{2,<:Real}` - A tuple of the (Input, Output) length scale factor. Default = (1,1)\\
`lB::NTuple{2,<:Real}` - A tuple of the (Input, Output) magnetic field scale factor. Default = (1,1)\\
`exp_mu0::NTuple{2,<:Real}` - A tuple of the (Input, Output) mu0 exponent (0, 1). Default = (0,0)
"""
function transform_cocos(cc_in::COCOS, cc_out::COCOS;
                         sigma_Ip::Union{NTuple{2,Int},Nothing} = nothing,
                         sigma_B0::Union{NTuple{2,Int},Nothing} = nothing,
                         ld::NTuple{2,<:Real} = (1,1),
                         lB::NTuple{2,<:Real} = (1,1),
                         exp_mu0::NTuple{2,<:Real} = (0,0))

    ld_eff = ld[2]/ld[1]
    lB_eff = lB[2]/lB[1]
    exp_mu0_eff = exp_mu0[2] - exp_mu0[1]

    sigma_RpZ_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ

    if sigma_Ip == nothing
        sigma_Ip_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    else
        sigma_Ip_eff = sigma_Ip[1]*sigma_Ip[2]
    end

    if sigma_B0 == nothing
        sigma_B0_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    else
        sigma_B0_eff = sigma_B0[1]*sigma_B0[2]
    end

    sigma_Bp_eff = cc_in.sigma_Bp * cc_out.sigma_Bp
    exp_Bp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    sigma_rhotp_eff = cc_in.sigma_rhotp * cc_out.sigma_rhotp

    mu0 = 4*pi*1e-7

    transforms = Dict()
    transforms["R"]        = ld_eff
    transforms["Z"]        = ld_eff
    transforms["P"]        = (lB_eff^2)/(mu0^exp_mu0_eff)
    transforms["PSI"]      = lB_eff * ld_eff^2 * sigma_Ip_eff * sigma_Bp_eff * ((2pi)^exp_Bp_eff) * ld_eff^2 * lB_eff
    transforms["ψ"]        = transforms["PSI"]
    transforms["TOR"]      = lB_eff * ld_eff^2 * sigma_B0_eff
    transforms["Φ"]        = transforms["TOR"]
    transforms["PPRIME"]   = (lB_eff/((ld_eff^2)*(mu0^exp_mu0_eff))) * sigma_Ip_eff * sigma_Bp_eff / ((2pi)^exp_Bp_eff)
    transforms["F_FPRIME"] = lB_eff * sigma_Ip_eff * sigma_Bp_eff / ((2pi)^exp_Bp_eff)
    transforms["B"]        = lB_eff * sigma_B0_eff
    transforms["F"]        = sigma_B0_eff * ld_eff * lB_eff
    transforms["I"]        = sigma_Ip_eff * ld_eff * lB_eff / (mu0^exp_mu0_eff)
    transforms["J"]        = sigma_Ip_eff * lB_eff/((mu0^exp_mu0_eff)*ld_eff)
    transforms["Q"]        = sigma_Ip_eff * sigma_B0_eff * sigma_rhotp_eff

    return transforms
end

function transform_cocos(cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS}; kwargs...)
    transform_cocos(cocos(cc_in), cocos(cc_out); kwargs...)
end

"""
    identify_cocos(sigma_B0, sigma_Ip, sigma_q, sigma_dpsi, clockwise_phi) -> List of possible COCOS IDs

Utility function to identify COCOS coordinate system. If multiple COCOS are possible, then all are returned.

Arguments:\\
`sigma_B0` - (+1,-1) sign of the toroidal magnetic field\\
`sigma_Ip` - (+1,-1) sign of toroidal plasma current\\
`sigma_q`  - (+1,-1) sign of the safety factor (q) within the plasma\\
`sigma_dpsi` -  +1 if psi is increasing, -1 if psi is decreasing\\
`clockwise_phi::Bool` - (optional) [True, False] if phi angle is defined clockwise or not. This is required to identify odd Vs even COCOS. Note that this cannot be determined from the output of a code.
"""
function identify_cocos(sigma_B0, sigma_Ip, sigma_q, sigma_dpsi,
                        clockwise_phi::Union{Bool,Nothing} = nothing)

    if clockwise_phi == nothing
        sigma_rpz = clockwise_phi
    elseif clockwise_phi
        sigma_rpz = -1
    else
        sigma_rpz = +1
    end

    # return both even and odd COCOS if clockwise_phi is not provided
    if sigma_rpz == nothing
        tmp  = identify_cocos(sigma_B0, sigma_Ip, sigma_q, sigma_dpsi, true)
        tmp2  = identify_cocos(sigma_B0, sigma_Ip, sigma_q, sigma_dpsi, false)
        return (tmp..., tmp2...)
    end

    sigma_Bp = sigma_dpsi / sigma_Ip
    sigma_rhotp = sigma_q / (sigma_Ip * sigma_B0)

    sigma2cocos = Dict(
        (+1, +1, +1) => 1,  # +Bp, +rpz, +rtp
        (+1, -1, +1) => 2,  # +Bp, -rpz, +rtp
        (-1, +1, -1) => 3,  # -Bp, +rpz, -rtp
        (-1, -1, -1) => 4,  # -Bp, -rpz, -rtp
        (+1, +1, -1) => 5,  # +Bp, +rpz, -rtp
        (+1, -1, -1) => 6,  # +Bp, -rpz, -rtp
        (-1, +1, +1) => 7,  # -Bp, +rpz, +rtp
        (-1, -1, +1) => 8)  # -Bp, -rpz, +rtp

    return (sigma2cocos[(sigma_Bp, sigma_rpz, sigma_rhotp)], sigma2cocos[(sigma_Bp, sigma_rpz, sigma_rhotp)] + 10)
end

"""
    identify_cocos(B0, Ip, q, psi, clockwise_phi, a) -> List of possible COCOS IDs

Utility function to identify COCOS coordinate system. If multiple COCOS are possible, then all are returned.

Arguments:\\
`B0` - toroidal magnetic field (with sign)\\
`Ip` - plasma current (with sign)\\
`q`  -  safety factor profile (with sign) as function of psi\\
`psi::AbstractVector` -  Vector of poloidal fluxs from the magnetic axis to the plasma boundary\\
`clockwise_phi::Bool` - (optional) [True, False] if phi angle is defined clockwise or not. This is required to identify odd Vs even COCOS. Note that this cannot be determined from the output of a code.\\
`a::AbstractVector` - (optional) flux surfaces minor radius as function of psi. This is required to identify 2*pi term in psi definition
"""
function identify_cocos(B0, Ip, q::AbstractVector, psi::AbstractVector,
                        clockwise_phi::Union{Bool,Nothing} = nothing,
                        a::Union{AbstractVector,Nothing} = nothing)

    sigma_B0 = sign(B0)
    sigma_Ip = sign(Ip)
    sigma_q = sign(q[1])
    sigma_dpsi = sign(psi[2]-psi[1])

    ccs = identify_cocos(sigma_B0, sigma_Ip, sigma_q, sigma_dpsi, clockwise_phi)

    # identify 2*pi term in psi definition based on q estimate
    if a != nothing
        index = argmin(abs.(q))
        if index == 1
            index += 1
        end

        q_estimate = abs.((pi * B0 * (a .- a[1]).^2) / (psi .- psi[1]))

        if abs(q_estimate[index] - q[index]) < abs(q_estimate[index] / (2 * pi) - q[index])
            eBp = 1
        else
            eBp = 0
        end

        return filter(x -> x > 10*eBp, ccs)
    else
        # return COCOS<10 as well as COCOS>10 if a is not provided
        return ccs
    end
end

# ----- GEQDSK Interface -----

"""
    check_cocos(g::GEQDSKFIle, cc::Union{Int,COCOS}) -> Bool

Checks if `GEQDSKFile` is consistant with the given cocos, `cc`
"""
function check_cocos(g::GEQDSKFile, cc::Union{Int,COCOS}; kwargs...)
    return check_cocos(g.bcentr, g.current, g.fpol, g.pprime, g.qpsi, g.psi, cocos(cc); kwargs...)
end

"""
    identify_cocos(g::GEQDSKFIle; clockwise_phi=nothing) -> List of possible COCOS IDs

Identifies possible `GEQDSKFile` COCOS. A unique identification requires setting the `clockwise_phi` keyword.
"""
function identify_cocos(g::GEQDSKFile; clockwise_phi = nothing)
    return filter(x -> x < 10, identify_cocos(g.bcentr, g.current, g.qpsi, g.psi, clockwise_phi, nothing))
end

"""
    cocos(g::GEQDSKFile; kwargs...) -> COCOS

Identifies and returns `GEQDSKFile` `COCOS`. A unique identification requires setting the `clockwise_phi` keyword.
"""
function cocos(g::GEQDSKFile; kwargs...)
    cc = identify_cocos(g; kwargs...)
    if length(cc) > 1
        @error "Unable to determine unique COCOS. Try providing clockwise_phi::Bool keyword. Possibilities are $cc"
    end
    return cocos(cc[1])
end

"""
    transform_cocos(g::GEQDSKFile, cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS}; kwargs...) -> GEQDSKFile

Transforms the given `GEQDSKFile` with `COCOS=cc_in`, and returns a `GEQDSKFile` with `COCOS=cc_out`.
"""
function transform_cocos(g::GEQDSKFile, cc_in::Union{Int,COCOS}, cc_out::Union{Int,COCOS}; kwargs...)
    T = transform_cocos(cc_in, cc_out; kwargs...)

    g_new = GEQDSKFile(g.file*" w/ cocos = $(cocos(cc_out).cocos)", g.nw, g.nh,
                       g.r*T["R"], g.z*T["Z"], g.rdim*T["R"], g.zdim*T["Z"],
                       g.rleft*T["R"], g.zmid*T["Z"], g.nbbbs, g.rbbbs*T["R"], g.zbbbs*T["Z"],
                       g.limitr, g.rlim*T["R"], g.zlim*T["Z"], g.rcentr*T["R"], g.bcentr*T["B"],
                       g.rmaxis*T["R"], g.zmaxis*T["Z"], g.simag*T["PSI"], g.sibry*T["PSI"],
                       g.psi*T["PSI"], g.current*T["I"], g.fpol*T["F"],
                       g.pres*T["P"], g.ffprim*T["F_FPRIME"], g.pprime*T["PPRIME"],
                       g.qpsi*T["Q"], g.psirz*T["PSI"])

    return g_new
end

