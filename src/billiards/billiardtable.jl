export Billiard, randominside
#######################################################################################
## Billiard Table
#######################################################################################
struct Billiard{T, D, O<:Tuple}
    obstacles::O
end

#pretty print:

_get_name(o::Obstacle) = :name ∈ fieldnames(typeof(o)) ? o.name : string(typeof(o))
function Base.show(io::IO, bd::Billiard{T,D,BT}) where {T, D, BT}
    s = "Billiard{$T} with $D obstacles:\n"
    for o in bd
        s*="  $(_get_name(o))\n"
    end
    print(io, s)
end



"""
    Billiard(obstacles...)
Construct a `Billiard` from given `obstacles` (tuple, vector, varargs).

For functions like [`boundarymap`](@ref),
it is expected (if possible) that the obstacles of the billiard are sorted,
such that the arc-coordinate `ξ` around the billiard is increasing counter-clockwise.

`ξ` is measured as:
* the distance from start point to end point in `Wall`s
* the arc length measured counterclockwise from the open face in `Semicircle`s
* the arc length measured counterclockwise from the rightmost point
  in `Circular`s
"""
function Billiard(bd::Union{AbstractVector, Tuple})

    T = eltype(bd[1])
    D = length(bd)
    # Assert that all elements of `bd` are of same type:
    for i in 2:D
        eltype(bd[i]) != T && throw(ArgumentError(
        "All obstacles of the billiard must have same type of
        numbers. Found $T and $(eltype(bd[i])) instead."
        ))
    end

    tup = (bd...,)
    return Billiard{T, D, typeof(tup)}(tup)
end

function Billiard(bd::Vararg{Obstacle})
    T = eltype(bd[1])
    tup = (bd...,)
    return Billiard(tup)
end


Base.getindex(bd::Billiard, i) = bd.obstacles[i]
Base.lastindex(bd::Billiard) = length(bd.obstacles)
Base.firstindex(bd::Billiard) = 1
# Iteration:
Base.iterate(bd::Billiard) = iterate(bd.obstacles)
Base.iterate(bd::Billiard, state) = iterate(bd.obstacles, state)
Base.length(::Billiard{T, L}) where {T, L} = L
Base.eltype(::Billiard{T}) where {T} = T

@inline isperiodic(bd) = count(x -> typeof(x) <: PeriodicWall, bd.obstacles) ≥ 2

#######################################################################################
## Distances
#######################################################################################
for f in (:distance, :distance_init)
    @eval $(f)(p::AbstractParticle, bd::Billiard) = $(f)(p.pos, bd.obstacles)
    @eval $(f)(pos::SV{T}, bd::Billiard) where {T} = $(f)(pos, bd.obstacles)
end

for f in (:distance, :distance_init)
    @eval begin
        function ($f)(p::SV{T}, bd::Tuple)::T where {T}
            dmin::T = T(Inf)
            for obst in bd
                d::T = distance(p, obst)
                d < dmin && (dmin = d)
            end
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
    bd::Union{Vector{<:Obstacle{T}}, Billiard{T}}) where {T<:AbstractFloat}

    xmin::T = ymin::T = T(Inf)
    xmax::T = ymax::T = T(-Inf)
    for obst ∈ bd
        xs::T, ys::T, xm::T, ym::T = cellsize(obst)
        xmin = xmin > xs ? xs : xmin
        ymin = ymin > ys ? ys : ymin
        xmax = xmax < xm ? xm : xmax
        ymax = ymax < ym ? ym : ymax
    end
    return xmin, ymin, xmax, ymax
end

"""
    randominside(bd::Billiard [, ω])
Return a particle with random allowed initial conditions inside the given
billiard. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω`.
"""
randominside(bd::Billiard) = Particle(_randominside(bd)...)
randominside(bd::Billiard{T}, ω) where {T} =
MagneticParticle(_randominside(bd)..., T(ω))



function _randominside(bd::Billiard{T}) where {T<:AbstractFloat}
    #1. position
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bd)

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    pos = SV{T}(xp, yp)

    dist = distance_init(pos, bd)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        pos = SV{T}(xp, yp)
        dist = distance_init(pos, bd)
    end

    #2. velocity
    φ = T(2π*rand())
    vel = SV{T}(sincos(φ)...)

    #3. current_cell (does nothing)
    cc = zero(SV{T})

    return pos, vel, cc
end
