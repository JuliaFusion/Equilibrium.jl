struct COCOS
    cocos::Int           # COCOS ID number
    exp_Bp::Int          # 0 or 1, depending if psi is already divided by 2pi or not, respectively
    sigma_Bp::Int        # +1 or -1, depending if psi is increasing or decreasing with Ip and B0 positive
    sigma_RpZ::Int       # +1 or -1, depending if (R,phi,Z) is right-handed or (R,Z,phi), respectively
    sigma_rhotp::Int     # +1 or -1, depending if (rho, theta, phi) is right-handed or (rho,phi,theta), repectively
    sign_q_pos::Int      # +1 or -1, depending if q is positive or negative with Ip and B0 positive
    sign_pprime_pos::Int # +1 or -1, depending if dp/dpsi is positive or negative with Ip and B0 positive
end

"""
Returns COCOS structure given the ID number
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
Returns True if GEQDSKFile is consistant with given COCOS
"""
function check(g::GEQDSKFile, cc::COCOS; verbose=false)

    valid = true
    qsign = sign(g.qpsi[end])
    if qsign*cc.sigma_rhotp*sign(g.current)*sign(g.bcentr) < 0
        verbose && @warn "sign(q[end]) ≠ sigma_rhotp*sign(Ip)*sign(B0)"
        valid = false
    end

    if all(sign.(g.fpol)*sign(g.bcentr) .< 0)
        verbose && @warn "Signs of F and B0 are not consistant"
        valid = false
    end

    if sign(g.sibry - g.simag)*cc.sigma_Bp*sign(g.current) < 0
        if g.sibry > g.simag
            verbose && @warn "psi should be decreasing with sign(Ip) = $(sign(g.current)) for COCOS = $(cc.cocos)"
        else
            verbose && @warn "psi should be increasing with sign(Ip) = $(sign(g.current)) for COCOS = $(cc.cocos)"
        end
        valid = false
    elseif sign(sum(g.pprime))*sign(g.current)*cc.sigma_Bp > 0
        verbose && @warn "sign(pprime) should be $(-sign(g.current)*cc.sigma_Bp)"
        valid = false
    end

    return valid
end

"""
Returns the native COCOS that an unmodified gEQDSK would obey, defined by sign(Bt) and sign(Ip)
In order for psi to increase from axis to edge and for q to be positive:
All use sigma_RpZ=+1 (phi is counterclockwise) and exp_Bp=0 (psi is flux/2.*pi)
We want
sign(psi_edge-psi_axis) = sign(Ip)*sigma_Bp > 0  (psi always increases in gEQDSK)
sign(q) = sign(Ip)*sign(Bt)*sigma_rhotp > 0      (q always positive in gEQDSK)
::
    ============================================
    Bt    Ip    sigma_Bp    sigma_rhotp    COCOS
    ============================================
    +1    +1       +1           +1           1
    +1    -1       -1           -1           3
    -1    +1       +1           -1           5
    -1    -1       -1           +1           7
"""
function cocos(B0,Ip)
    cc = cocos(1)
    sign_Bt = Int(cc.sigma_RpZ*sign(B0))
    sign_Ip = Int(cc.sigma_RpZ*sign(Ip))

    g_cocos = Dict((+1, +1) => 1, # +Bt, +Ip
                   (+1, -1) => 3, # +Bt, -Ip
                   (-1, +1) => 5, # -Bt, +Ip
                   (-1, -1) => 7, # -Bt, -Ip
                   (+1,  0) => 1, # +Bt, No current
                   (-1,  0) => 3) # -Bt, No current

    return g_cocos[(sign_Bt,sign_Ip)]
end

function cocos(g::GEQDSKFile)
    return cocos(g.bcentr,g.current)
end

function transforms(cc_in::COCOS, cc_out::COCOS)
    sigma_Ip_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    sigma_B0_eff = cc_in.sigma_RpZ * cc_out.sigma_RpZ
    sigma_Bp_eff = cc_in.sigma_Bp * cc_out.sigma_Bp
    exp_Bp_eff = cc_out.exp_Bp - cc_in.exp_Bp
    sigma_rhotp_eff = cc_in.sigma_rhotp * cc_out.sigma_rhotp

    transforms = Dict()
    transforms["1/PSI"] = sigma_Ip_eff * sigma_Bp_eff / ((2*pi)^exp_Bp_eff)
    transforms["invPSI"] = transforms["1/PSI"]
    transforms["dPSI"] = transforms["1/PSI"]
    transforms["F_FPRIME"] = transforms["dPSI"]
    transforms["PPRIME"] = transforms["dPSI"]
    transforms["PSI"] = sigma_Ip_eff * sigma_Bp_eff * ((2pi)^exp_Bp_eff)
    transforms["Q"] = sigma_Ip_eff * sigma_B0_eff * sigma_rhotp_eff
    transforms["TOR"] = sigma_B0_eff
    transforms["BT"] = transforms["TOR"]
    transforms["IP"] = transforms["TOR"]
    transforms["F"] = transforms["TOR"]
    transforms["POL"] = sigma_B0_eff * sigma_rhotp_eff
    transforms["BP"] = transforms["POL"]

    return transforms
end


function transforms(cc_in::Int, cc_out::Int)
    transforms(cocos(cc_in), cocos(cc_out))
end
