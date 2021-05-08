struct Boundary{T}
    points::Vector{NTuple{2,T}}

    function Boundary(points::Vector{NTuple{2,T}}) where T
        if first(points) != last(points)
            push!(points,first(points))
        end
        return new{T}(points)
    end
end

const PlasmaBoundary = Boundary
const FluxSurface = Boundary
const Wall = Boundary

function Base.show(io::IO, b::Boundary)
    print(io,"Boundary: npoints = $(length(b.points))")
end

Base.broadcastable(b::Boundary) = (b,)

function Boundary(x::Vector,y::Vector)
    Boundary(collect(zip(x,y)))
end

function Base.getproperty(b::Boundary,s::Symbol)
    if s == :points
        return getfield(b,s)
    end
    if s in (:r,:x)
        v = getfield(b,:points)
        return getindex.(v,1)
    end
    if s in (:z,:y)
        v = getfield(b,:points)
        return getindex.(v,2)
    end
    error("type $(typeof(b)) has no field $s")
end

in_boundary(b::Boundary, p) = PolygonOps.inpolygon(p,b.points,in=true,on=true,out=false)
in_boundary(b::Boundary,x,y) = in_boundary(b,(x,y))

const in_vessel = in_boundary
const in_plasma = in_boundary

function boundary(M::AbstractEquilibrium, psi, dx=0.01,dy=0.01)
    xlims, ylims = limits(M)
    x = range(xlims...,step=dx)
    y = range(ylims...,step=dy)
    Psi = [M(xx,yy) for xx in x, yy in y]
    cl = Contour.contour(x,y,Psi,psi)
    l = Contour.lines(cl)[1]
    xc, yc = Contour.coordinates(l)
    return Boundary(xc,yc)
end

function circumference(b::Boundary)
    p = b.points
    return sum(norm(p[i+1] .- p[i]) for i=1:(length(p)-1))
end

function average(b::Boundary, F)
    p = b.points
    return sum(norm(p[i+1] .- p[i])*(F(p[i][1],p[i][2])) for i=1:(length(p)-1))/circumference(b)
end

area(s::Boundary) = PolygonOps.area(s.points)

function surface_area(b::Boundary)
    p = b.points
    return 2pi*sum(norm(p[i+1] .- p[i])*(p[i+1][1] + p[i][1])/2 for i=1:(length(p)-1))
end

function surface_area_average(b::Boundary, F)
    A = surface_area(b)
    p = b.points
    return (2pi/A)*sum(norm(p[i+1] .- p[i])*((p[i+1][1] + p[i][1])/2)*
                   (F(p[i+1][1],p[i+1][2]) + F(p[i][1],p[i][2]))/2 for i=1:(length(p)-1))
end

function surface_area_average(b::Boundary, F::Vector)
    A = surface_area(b)
    p = b.points
    return (2pi/A)*sum(norm(p[i+1] .- p[i])*((p[i+1][1] + p[i][1])/2)*
                       (F[i+1] + F[i])/2 for i=1:(length(p)-1))
end

function area_average(b::Boundary, F; dx=0.01, dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    A = zero(dx)
    for xx in x, yy in y
        if in_boundary(b,(xx,yy))
            A = A + F(xx,yy)*dx*dy
        end
    end
    return A/area(b)
end

function area_average(b::Boundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    dx = step(x)
    dy = step(y)
    A = zero(eltype(x))
    for ix in eachindex(x), iy in eachindex(y)
        if in_boundary(b, (x[ix],y[iy]))
            A = A + F[ix,iy]*dx*dy
        end
    end
    return A/area(b)
end

function volume(b::Boundary; dx=0.01,dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    V = zero(eltype(x))
    for xx in x, yy in y
        if in_boundary(b, (xx,yy))
            V = V + xx*dx*dy
        end
    end
    return 2pi*V
end

function volume_average(b::Boundary, F; dx=0.01, dy=0.01)
    x = range(extrema(b.r)...,step=dx) .+ dx/2
    y = range(extrema(b.z)...,step=dy) .+ dy/2

    A = zero(dx)
    for xx in x, yy in y
        if in_boundary(b, (xx,yy))
            A = A + F(xx,yy)*xx*dx*dy
        end
    end
    return 2pi*A/volume(b)
end

function volume_average(b::Boundary, F::AbstractMatrix, x::AbstractRange, y::AbstractRange)
    dx = step(x)
    dy = step(y)
    A = zero(eltype(x))
    for ix in eachindex(x), iy in eachindex(y)
        if in_boundary(b, (x[ix],y[iy]))
            A = A + F[ix,iy]*x[ix]*dx*dy
        end
    end
    return 2pi*A/volume(b)
end
