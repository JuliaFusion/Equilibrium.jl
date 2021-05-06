struct PlasmaBoundary{T}
    points::Vector{NTuple{2,T}}
end

function Base.show(io::IO, l::PlasmaBoundary)
    print(io,"PlasmaBoundary: npoints = $(length(l.points))")
end

Base.broadcastable(l::PlasmaBoundary) = (l,)

function PlasmaBoundary(points)
    if first(points) != last(points)
        append!(points,first(points))
    end
    return PlasmaBoundary(points)
end

function PlasmaBoundary(x::Vector,y::Vector)
    PlasmaBoundary(collect(zip(x,y)))
end

function Base.getproperty(sh::PlasmaBoundary{T},s::Symbol) where T<:Real
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

in_plasma(s::PlasmaBoundary, p) = inpolygon(p,s.points,in=true,on=true,out=false)
in_plasma(s::PlasmaBoundary,x,y) = in_plasma(s,(x,y))

function circumference(s::PlasmaBoundary)
    return sum(norm(s.points[i+1] .- s.points[i]) for i=1:(length(s.points)-1))
end

area(s::PlasmaBoundary) = PolygonOps.area(s.points)

function area_average(s::PlasmaBoundary, F; dx=0.01, dy=0.01)
    x = range(extrema(s.r)...,step=dx) .+ dx/2
    y = range(extrema(s.z)...,step=dy) .+ dy/2

    A = zero(dx)
    for xx in x, yy in y
        if inpolygon((xx,yy),s.points)
            A = A + F(xx,yy)*dx*dy
        end
    end
    return A
end

function area_average(s::PlasmaBoundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    dx = step(x)
    dy = step(x)
    A = zero(rmin)
    for ix in eachindex(x), iy in eachindex(y)
        if inpolygon((x[ix],y[iy]),s.points)
            A = A + F[ix,iy]*dx*dy
        end
    end
    return A
end

function volume(s::PlasmaBoundary;dx=0.01,dy=0.01)
    x = range(extrema(s.r)...,step=dx) .+ dx/2
    y = range(extrema(s.z)...,step=dy) .+ dy/2

    V = zero(rmin)
    for xx in x, yy in y
        if inpolygon((xx,yy),s.points)
            V = V + xx*dx*dy
        end
    end
    return 2pi*V
end

function volume_average(s::PlasmaBoundary, F; dx=0.01, dy=0.01)
    x = range(extrema(s.r)...,step=dx) .+ dx/2
    y = range(extrema(s.z)...,step=dy) .+ dy/2

    A = zero(dx)
    for xx in x, yy in y
        if inpolygon((xx,yy),s.points)
            A = A + F(xx,yy)*xx*dx*dy
        end
    end
    return 2pi*A
end

function volume_average(s::PlasmaBoundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    dx = step(x)
    dy = step(x)
    A = zero(rmin)
    for ix in eachindex(x), iy in eachindex(y)
        if inpolygon((x[ix],y[iy]),s.points)
            A = A + F[ix,iy]*x[ix]*dx*dy
        end
    end
    return 2pi*A
end
