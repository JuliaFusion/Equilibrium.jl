using EFIT

function load(g::GEQDSKFile)
    r = linspace(extrema(g.r)...,length(g.r))
    z = linspace(extrema(g.z)...,length(g.z))
    psi = linspace(extrema(g.psi)...,length(g.psi))
    M = AxisymmetricEquilibrium(r, z, psi, g.psirz, g.fpol, g.pres, g.qpsi, g.fpol*0, (g.rmaxis,g.zmaxis),g.sibry);
    return M
end

function load_geqdsk(gfile)
    return load(readg(gfile))
end

function load_limiter(g::GEQDSKFile)
    lim = Limiter()
    for i =1:length(g.rlim)
        push!(lim.vertices,(g.rlim[i],g.zlim[i]))
    end
    push!(lim.vertices, (g.rlim[1],g.zlim[1]))

    return lim
end

function load_geqdsk_limiter(gfile)
    return load_limiter(readg(gfile))
end
