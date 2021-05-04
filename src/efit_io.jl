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

function Limiter(g::GEQDSKFile)
    lim = Limiter()
    for i =1:length(g.rlim)
        push!(lim.vertices,(g.rlim[i],g.zlim[i]))
    end
    push!(lim.vertices, (g.rlim[1],g.zlim[1]))

    return lim
end

function read_geqdsk(gfile; kwargs...)
    g = readg(gfile)
    return AxisymmetricEquilibrium(g; kwargs...), Limiter(g)
end
