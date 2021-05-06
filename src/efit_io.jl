using EFIT

function AxisymmetricEquilibrium(g::GEQDSKFile, cc::Union{Int,COCOS})
    M = AxisymmetricEquilibrium(cocos(cc), g.r, g.z, g.psi, g.psirz, g.fpol, g.pres,
                                g.qpsi, g.fpol*0, (g.rmaxis,g.zmaxis));
    return M
end

function AxisymmetricEquilibrium(g::GEQDSKFile; clockwise_phi::Bool)
    cc = cocos(g,clockwise_phi=clockwise_phi)
    M = AxisymmetricEquilibrium(g, cc)
    return M
end

function Wall(g::GEQDSKFile)
    return Wall(collect(zip(g.rlim,g.zlim)))
end

function read_geqdsk(gfile; kwargs...)
    g = readg(gfile)
    return AxisymmetricEquilibrium(g; kwargs...), Wall(g)
end
