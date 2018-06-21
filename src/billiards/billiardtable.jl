export Billiard, randominside
import Base:iterate
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

If you want to use the [`boundarymap`](@ref) function, then it is expected to
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
iterate(bt::Billiard) = iterate(bt.obstacles)
iterate(bt::Billiard, state) = iterate(bt.obstacles, state)

eltype(bt::Billiard{T}) where {T} = T

isperiodic(bt) = any(x -> typeof(x) <: PeriodicWall, bt.obstacles)

# total arclength
totallength(bt::Billiard) = sum(totallength(x) for x in bt.obstacles)

#######################################################################################
## Distances
#######################################################################################
for f in (:distance, :distance_init)
    @eval $(f)(p::AbstractParticle, bt::Billiard) = $(f)(p.pos, bt.obstacles)
    @eval $(f)(pos::SV{T}, bt::Billiard) where {T} = $(f)(pos, bt.obstacles)
end

for f in (:distance, :distance_init)
    @eval begin
        function ($f)(p::SV{T}, bt::Tuple)::T where {T}
            dmin::T = T(Inf)
            for obst in bt
                d::T = distance(p, obst)
                d < dmin && (dmin = d)
            end#obstacle loop
            return dmin
        end
    end
end

#######################################################################################
## total arclength
#######################################################################################
function totallength(bt::Billiard)
    #for some reason, this is faster than @inline totallength(bt) = ...
    return unrolled_reduce(+,0.0, unrolled_map(x->totallength(x),bt.obstacles))
end

function real_coordinates(ξ, sφ, bt::Billiard{T}; return_obst::Bool = false) where T
    abs(sφ) > 1 && throw(DomainError())#"|sin φ| must not be larger than 1"))
    lower = zero(T)
    upper = lower
    for (i, obst) ∈ enumerate(bt)
        #println("testing $(obst.name)")
        upper = lower + totallength(obst)
        #println("\tbounds: $lower:$upper")
        if ξ <= upper
            ret = real_coordinates(ξ - lower, sφ, obst)
            return return_obst ? (ret..., i) : ret
        end
        lower = upper
        #println("\tNEXT!")
    end
    #println("was: $ξ\tmax: $upper")
    throw(DomainError())#"ξ is too large for this billiard!"))
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
randominside(bt::Billiard) = Particle(_randominside(bt)...)
randominside(bt::Billiard{T}, ω) where {T} =
MagneticParticle(_randominside(bt)..., T(ω))



function _randominside(bt::Billiard{T}) where {T<:AbstractFloat}
    #1. position
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

    #2. velocity
    φ = T(2π*rand())
    vel = SV{T}(sin(φ), cos(φ)) #TODO:Change to sincos for julia 0.7

    #3. current_cell (does nothing)
    cc = zero(SV{T})

    return pos, vel, cc
end
