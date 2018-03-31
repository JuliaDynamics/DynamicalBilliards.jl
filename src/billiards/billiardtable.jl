export Billiard, randominside
import Base:start, next, done
#######################################################################################
## Billiard Table
#######################################################################################
struct Billiard{T, D, O<:Tuple}
    obstacles::O
end

#pretty print:
function Base.show(io::IO, bt::Billiard{T,D,BT}) where {T, D, BT}
    s = "Billiard{$T} with $D obstacles:\n"
    for o in bt
        s*="  $(o.name)\n"
    end
    print(io, s)
end



"""
    Billiard(obstacles...)
Construct a `Billiard` from given `obstacles` (tuple, vector, varargs).

If you want to use the [`poincaresection`](@ref) function, then it is expected to
provide the obstacles of the billiard in sorted order, such that the boundary
coordinate (measured using [`arclength`](@ref))
around the billiard is increasing counter-clockwise.

The boundary coordinate is measured as:
* the distance from start point to end point in `Wall`s
* the arc length measured counterclockwise from the open face in `Semicircle`s
* the arc length measured counterclockwise from the rightmost point in `Circular`s
"""
function Billiard(bt::Union{AbstractVector, Tuple})

    T = eltype(bt[1])
    D = length(bt)
    # Assert that all elements of `bt` are of same type:
    for i in 2:D
        eltype(bt[i]) != T && throw(ArgumentError(
        "All obstacles of the billiard must have same type of
        numbers. Found $T and $(eltype(bt[i])) instead."
        ))
    end

    tup = (bt...,)
    return Billiard{T, D, typeof(tup)}(tup)
end

function Billiard(bt::Vararg{Obstacle})
    T = eltype(bt[1])
    tup = (bt...,)
    return Billiard(tup)
end


getindex(bt::Billiard, i) = bt.obstacles[i]
# Iteration:
start(bt::Billiard) = start(bt.obstacles)
next(bt::Billiard, state) = next(bt.obstacles, state)
done(bt::Billiard, state) = done(bt.obstacles, state)

eltype(bt::Billiard{T}) where {T} = T


isperiodic(bt) = Unrolled.unrolled_any(x -> typeof(x) <: PeriodicWall, bt)


#######################################################################################
## Distances
#######################################################################################
for f in (:distance, :distance_init)
    @eval $(f)(p::AbstractParticle, bt::Billiard) = $(f)(p.pos, bt.obstacles)
end

for f in (:distance, :distance_init)
    @eval begin
        @unroll function ($f)(p::SV{T}, bt::Tuple)::T where {T}
            dmin::T = T(Inf)
            @unroll for obst in bt
                d::T = distance(p, obst)
                d < dmin && (dmin = d)
            end#obstacle loop
            return dmin
        end
    end
end

function distance(p::AbstractParticle{T}, bt::Billiard{T})::T where {T}
    d = T(Inf)
    for obst ∈ bt
        di = distance(p, obst)
        di < d && (d = di)
    end
    return d
end

function distance_init(pos::SVector, bt::Billiard{T})::T where {T}
    d = T(Inf)
    for obst ∈ bt
        di = distance_init(pos, obst)
        di < d && (d = di)
    end
    return d
end


#######################################################################################
## randominside
#######################################################################################
function cellsize(
    bt::Union{Vector{<:Obstacle{T}}, Billiard{T}}) where {T<:AbstractFloat}

    xmin::T = ymin::T = T(Inf)
    xmax::T = ymax::T = T(-Inf)
    for obst ∈ bt
        xs::T, ys::T, xm::T, ym::T = cellsize(obst)
        xmin = xmin > xs ? xs : xmin
        ymin = ymin > ys ? ys : ymin
        xmax = xmax < xm ? xm : xmax
        ymax = ymax < ym ? ym : ymax
    end
    return xmin, ymin, xmax, ymax
end

"""
    randominside(bt::Billiard [, ω])
Return a particle with random allowed initial conditions inside the given
billiard. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω`.
"""
randominside(bt::Billiard{T}) where {T} =
    Particle(_randominside(bt)..., T(2π*rand()))
randominside(bt::Billiard{T}, ω) where {T} =
    MagneticParticle(_randominside(bt)..., T(2π*rand()), T(ω))

function _randominside(bt::Billiard{T}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    pos = SV{T}(xp, yp)

    dist = distance_init(pos, bt)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        pos = SV{T}(xp, yp)
        dist = distance_init(pos, bt)
    end

    return pos[1], pos[2]
end
