struct Wall{T,N}
    points::Vector{NTuple{N,T}}
end

const Wall2D{T} = Wall{T,2} where T

function Wall(points)
    if first(points) != last(points)
        append!(points,first(points))
    end
    return Wall(points)
end

function Base.show(io::IO, l::Wall)
    print(io,"Wall: npoints = $(length(l.points))")
end

Base.broadcastable(l::Wall) = (l,)

function Base.getproperty(l::Wall2D{T},s::Symbol) where T<:Real
    if s == :points
        return getfield(l,s)
    end
    if s == :r
        v = getfield(l,:points)
        return getindex.(v,1)
    end
    if s == :z
        v = getfield(l,:points)
        return getindex.(v,2)
    end
    error("type $(typeof(l)) has no field $s")
end

in_vessel(wall,p) = inpolygon(p, wall.points, in=true,on=false,out=false)
in_vessel(wall,r,z) = in_vessel(wall,(r,z))
