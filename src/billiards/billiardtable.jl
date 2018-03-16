export BilliardTable, randominside
####################################################
## Billiard Table
####################################################
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






####################################################
#%% randominside
"""
    randominside(bt::Vector{<:Obstacle{T}}[, ω])
Return a particle with correct (allowed) initial conditions inside the given
billiard table defined by the vector `bt`. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω` (unless
`ω` is 0). Else, it is `Particle`.
"""
function randominside(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = T(f*2π)

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = Particle([xp, yp, φ0])

    dist = distance_init(p, bt)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        p.pos = SVector{2,T}(xp, yp)
        dist = distance_init(p, bt)
    end

    return p
end

function randominside(ω::Real, bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    ω = convert(T, ω)
    if ω == 0
        return randominside(bt)
    end

    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = T(f*2π)

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = MagneticParticle([xp, yp, φ0], ω)

    dist = distance_init(p, bt)
    while dist <= sqrt(eps(T))

        xp = rand()*(xmax-xmin) + xmin
        yp = rand()*(ymax-ymin) + ymin
        p.pos = [xp, yp]
        dist = distance_init(p, bt)
    end

    return p
end
randominside(bt::Vector{<:Obstacle{T}}, ω::Real) where {T} = randominside(ω, bt)
