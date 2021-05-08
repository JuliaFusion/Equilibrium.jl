using EFIT

function efit(g::GEQDSKFile, cc::Union{Int,COCOS})
    M = efit(cocos(cc), g.r, g.z, g.psi, g.psirz, g.fpol, g.pres,
             g.qpsi, g.fpol*0, (g.rmaxis,g.zmaxis),
             Int(sign(g.bcentr)*sign(g.current)));
    return M
end

function efit(g::GEQDSKFile; clockwise_phi::Bool)
    cc = cocos(g,clockwise_phi=clockwise_phi)
    M = efit(g, cc)
    return M
end

function Wall(g::GEQDSKFile)
    return Wall(collect(zip(g.rlim,g.zlim)))
end

function read_geqdsk(gfile; kwargs...)
    g = readg(gfile)
    return efit(g; kwargs...), Wall(g)
end
