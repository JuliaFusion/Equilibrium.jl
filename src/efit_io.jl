using EFIT

function AxisymmetricEquilibrium(g::GEQDSKFile)
    M = AxisymmetricEquilibrium(g.r, g.z, g.psi, g.psirz, g.fpol, g.pres, g.qpsi, g.fpol*0, (g.rmaxis,g.zmaxis),g.sibry);
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

function read_geqdsk(gfile)
    g = readg(gfile)
    return AxisymmetricEquilibrium(g), Limiter(g)
end
