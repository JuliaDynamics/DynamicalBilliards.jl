export escapetime, meancollisiontime
export escapetime!, meancollisiontime!
export visited_obstacles, visited_obstacles!

"""
    escapetime([p,] bd, t; warning = false)
Calculate the escape time of a particle `p` in the billiard `bd`, which
is the time until colliding with any "door" in `bd`.
As a "door" is considered any [`FiniteWall`](@ref) with
field `isdoor = true`.

If the particle evolves for more than `t` (integer or float) without colliding with the
`Door` (i.e. escaping) the returned result is `Inf`.

A warning can be thrown if the result is `Inf`. Enable this using the keyword
`warning = true`.

If a particle is not given, a random one is picked through [`randominside`](@ref).
See [`parallelize`](@ref) for a parallelized version.
"""
escapetime(p, bd, t; warning = false) =
    escapetime!(copy(p), bd, t; warning = warning)
escapetime(bd::Billiard, t; kwargs...) =
    escapetime(randominside(bd), bd, t; kwargs...)

function escapetime!(p::AbstractParticle{T}, bd::Billiard{T}, t, raysplitters = nothing;
        warning::Bool=false)::T where {T<:AbstractFloat}

    ei = escapeind(bd)
    if length(ei) == 0
        error("""The billiard does not have any "doors"!""")
    end

    totalt = zero(T); count = zero(t); t_to_write = zero(T)
    isray = !isa(raysplitters, Nothing)
    isray && acceptable_raysplitter(raysplitteers, bd)
    raysidx = raysplit_indices(bd, raysplitters)

    while count < t

        i, tmin, pos, vel = bounce!(p, bd, raysidx, raysplitters)
        t_to_write += tmin

        if isperiodic(i, bd)
            continue
        else
            totalt += t_to_write
            i ∈ ei &&  break # the collision happens with a Door!

            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
    if count ≥ t
        warning && warn("Particle did not escape within max-iter window.")
        return T(Inf)
    end
    return totalt
end

function escapeind(bd)
    j = Int[]
    for (i, obst) in enumerate(bd)
        if typeof(obst) <: FiniteWall && obst.isdoor == true
            push!(j, i)
        end
    end
    return j
end


"""
    meancollisiontime([p,] bd, t) → κ
Compute the mean collision time `κ` of the particle `p` in the billiard `bd` by
evolving for total amount `t` (either float for time or integer for collision number).

Collision times are counted only between obstacles that are *not*
[`PeriodicWall`](@ref).

If a particle is not given, a random one is picked through [`randominside`](@ref).
See [`parallelize`](@ref) for a parallelized version.
"""
meancollisiontime(p, bd, t) = meancollisiontime!(copy(p), bd, t)
meancollisiontime(bd::Billiard, t) = meancollisiontime(randominside(bd), bd, t)

function meancollisiontime!(p::AbstractParticle{T}, bd::Billiard{T}, t)::T where {T}

    ispinned(p, bd) && return Inf
    ismagnetic = typeof(p) <: MagneticParticle
    ismagnetic && (absω = abs(p.omega))
    tmin, i = next_collision(p, bd)
    tmin == Inf && return Inf

    κ = zero(T); count = zero(t); t_to_write = zero(T); colcount = 0

    while count < t
        i, tmin, pos, vel = bounce!(p, bd)

        t_to_write += tmin

        if isperiodic(bd) && i ∈ bd.peridx
            continue
        else
            κ += t_to_write
            colcount += 1
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end
    return κ/colcount
end

"""
    visited_obstacles(p, args...)
Same as [`visited_obstacles!`](@ref) but copies the particle instead.
"""
visited_obstacles(p::AbstractParticle, args...) =
    visited_obstacles!(copy(p), args...)
visited_obstacles(bd::Billiard, args...; kwargs...) =
    visited_obstacles!(randominside(bd), bd, args...; kwargs...)


"""
    visited_obstacles!([p::AbstractParticle,] bd::Billiard, t)
Evolve the given particle `p` inside the billiard `bd` exactly like
[`evolve!`](@ref). However return only:

* `ts::Vector{T}` : Vector of time points of when each collision occured.
* `obst::Vector{Int}` : Vector of obstacle indices in `bd` that the particle
  collided with at the time points in `ts`.

The first entries are `0.0` and `0`.
Similarly with [`evolve!`](@ref) the function does not
record collisions with periodic walls.

Currently does not support raysplitting. Returns empty arrays
for pinned particles.
"""
function visited_obstacles!(
    p::AbstractParticle{T}, bd::Billiard{T}, t) where {T<:AbstractFloat}

    if t ≤ 0
        throw(ArgumentError("cannot evolve backwards in time."))
    end
    if ispinned(p, bd)
        return T[], Int[]
    end

    ts = T[0.0]; obst = Int[0]

    count = zero(t); t_to_write = zero(T)
    if typeof(t) == Int
        for zzz in (ts, obst)
            sizehint!(zzz, t)
        end
    end

    while count < t
        i, tmin, pos, vel = bounce!(p, bd)
        t_to_write += tmin
        if isperiodic(i, bd)
            continue
        else
            count += increment_counter(t, tmin)
            push!(ts, t_to_write + ts[end]); push!(obst, i)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
    return ts, obst
end
