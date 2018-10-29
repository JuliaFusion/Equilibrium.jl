struct Limiter{T}
    vertices::Vector{NTuple{2,T}}
end

Limiter(T=Float64) = Limiter(NTuple{2,T}[])

function Base.show(io::IO, l::Limiter)
    print(io,"Limiter: npoints = $(length(l.vertices))")
end

function Base.getproperty(l::Limiter{T},s::Symbol) where T<:Real
    if s == :vertices
        return getfield(l,s)
    end
    if s == :r
        v = getfield(l,:vertices)
        return getindex.(v,1)
    end
    if s == :z
        v = getfield(l,:vertices)
        return getindex.(v,2)
    end
    error("type $(typeof(l)) has no field $s")
end

function is_left(p0,p1,p2)
    return ((p1[1] - p0[1]) * (p2[2] - p0[2]) - (p2[1] -  p0[1]) * (p1[2] - p0[2]))
end

function in_vessel(lim,p)

    v = lim.vertices
    nv = length(v)
    wn = 0
    @inbounds for i = 1:nv-1
        if v[i][2] <= p[2]
            if v[i+1][2] > p[2]
                if is_left(v[i],v[i+1],p) > 0.0
                    wn = wn +1
                end
            end
        else
            if v[i+1][2] <= p[2]
                if is_left(v[i],v[i+1],p) < 0.0
                    wn = wn - 1
                end
            end
        end
    end
    if wn == 0
        return false
    else
        return true
    end
end

function in_vessel(lim,r,z)
    in_vessel(lim,[r,z])
end
