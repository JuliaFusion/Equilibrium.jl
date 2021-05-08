using EFIT

function NumericalEquilibrium(g::GEQDSKFile, cc::Union{Int,COCOS})
    M = NumericalEquilibrium(cocos(cc), g.r, g.z, g.psi, g.psirz, g.fpol, g.pres,
                             g.qpsi, g.fpol*0, (g.rmaxis,g.zmaxis),
                             Int(sign(g.bcentr)*sign(g.current)));
    return M
end

function NumericalEquilibrium(g::GEQDSKFile; clockwise_phi::Bool)
    cc = cocos(g,clockwise_phi=clockwise_phi)
    M = NumericalEquilibrium(g, cc)
    return M
end

function Wall(g::GEQDSKFile)
    return Wall(collect(zip(g.rlim,g.zlim)))
end

function read_geqdsk(gfile; kwargs...)
    g = readg(gfile)
    return NumericalEquilibrium(g; kwargs...), Wall(g)
end
