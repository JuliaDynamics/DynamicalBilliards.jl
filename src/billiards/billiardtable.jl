export BilliardTable, randominside
#######################################################################################
## Billiard Table
#######################################################################################
immutable BilliardTable{T, BT<:Tuple}
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



#######################################################################################
## next_collision
#######################################################################################
"""
    next_collision(p, bt) -> (tmin, index)
Return the minimum collision time out of all `collisiontime(p, obst)` for `obst ∈ bt`,
as well as the `index` of the corresponding obstacle.
"""
function next_collision(
    p::AbstractParticle{T}, bt::Vector{<:Obstacle{T}})::Tuple{T,Int} where {T}
    tmin::T = T(Inf)
    ind::Int = 0
    for i in eachindex(bt)
        tcol::T = collisiontime(p, bt[i])
        # Set minimum time:
        if tcol < tmin
            tmin = tcol
            ind = i
        end
    end#obstacle loop
    return tmin, ind
end

function next_collision(
    p::AbstractParticle{T}, bt::Tuple)::Tuple{T,Int} where {T}
    findmin(map(x -> collisiontime(p, x), bt))
end

##

#######################################################################################
## randominside
#######################################################################################
function cellsize(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}

    xmin::T = ymin::T = T(Inf)
    xmax::T = ymax::T = T(-Inf)
    for i in eachindex(bt)
        xs::T, ys::T, xm::T, ym::T = cellsize(bt[i])
        xmin = xmin > xs ? xs : xmin
        ymin = ymin > ys ? ys : ymin
        xmax = xmax < xm ? xm : xmax
        ymax = ymax < ym ? ym : ymax
    end
    return xmin, ymin, xmax, ymax
end

"""
    randominside(bt::Vector{<:Obstacle{T}}[, ω])
Return a particle with correct (allowed) initial conditions inside the given
billiard table defined by the vector `bt`. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω` (unless
`ω` is 0). Else, it is `Particle`.
"""
randominside(bt) = Particle(_randominside(bt)...)
randominside(bt, ω) = MagneticParticle(_randominside(bt)..., ω)

function _randominside(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = T(f*2π)

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

    return pos[1], pos[2], φ0
end
