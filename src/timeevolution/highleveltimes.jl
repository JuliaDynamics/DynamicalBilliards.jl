export escapetime, meancollisiontime
export escapetime!, meancollisiontime!

#######################################################################################
## Escape times
#######################################################################################
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

function escapetime!(p::AbstractParticle{T}, bd::Billiard{T}, t;
        warning::Bool=false)::T where {T<:AbstractFloat}

    ei = escapeind(bd)
    if length(ei) == 0
        error("""The billiard does not have any "doors"!""")
    end

    totalt = zero(T); count = zero(t); t_to_write = zero(T)

    while count < t

        i, tmin, pos, vel = bounce!(p, bd)
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

#######################################################################################
## Mean Collision Time
#######################################################################################
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
