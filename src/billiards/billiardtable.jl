export BilliardTable, randominside, isperiodic
#######################################################################################
## Billiard Table
#######################################################################################
struct BilliardTable{T, BT<:Tuple}
    bt::BT
end

function BilliardTable(bt)
    T = eltype(bt[1])
    if typeof(bt) <: Tuple
        return BilliardTable{T, typeof(bt)}(bt)
    else
        tup = (bt...)
        return BilliardTable{T, typeof(tup)}(tup)
    end
end

# Need to define iteration in billiard table (for obst in bt...)

getobstacle(bt::BilliardTable{T,S}, ::Val{N}) where {T,S,N} = bt.bt[N]

getindex(bt::BilliardTable, i) = bt.bt[i]

isperiodic(bt) = Unrolled.unrolled_any(x -> typeof(x) <: PeriodicWall, bt)


#######################################################################################
## Distances
#######################################################################################
distance(p::AbstractParticle, bt::BilliardTable) = distance(p.pos, bt.bt)
distance_init(p::AbstractParticle, bt::BilliardTable) = distance_init(p.pos, bt.bt)

@unroll function distance_init(p::SV{T}, bt::Tuple)::T where {T}
    dmin::T = T(Inf)
    @unroll for obst in bt
        d::T = distance_init(p, obst)
        d < dmin && (dmin = d)
    end#obstacle loop
    return dmin
end
@unroll function distance(p::SV{T}, bt::Tuple)::T where {T}
    dmin::T = T(Inf)
    @unroll for obst in bt
        d::T = distance(p, obst)
        d < dmin && (dmin = d)
    end#obstacle loop
    return dmin
end


#######################################################################################
## randominside
#######################################################################################
function cellsize(
    bt::Union{Vector{<:Obstacle{T}}, BilliardTable{T}}) where {T<:AbstractFloat}

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
    randominside(bt::Vector{<:Obstacle{T}} [, ω])
Return a particle with allowed initial conditions inside the given
billiard table. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω`.
"""
randominside(bt::Vector{<:Obstacle{T}}) where {T} =
    Particle(_randominside(bt)..., T(2π*rand()))
randominside(bt::Vector{<:Obstacle{T}}, ω) where {T} =
    MagneticParticle(_randominside(bt)..., T(2π*rand()), T(ω))

function _randominside(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
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
