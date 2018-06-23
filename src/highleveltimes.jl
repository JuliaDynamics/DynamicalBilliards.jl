export escapetime, meancollisiontime
export escapetime!, meancollisiontime!

#######################################################################################
## Escape times
#######################################################################################
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
    escapetime(p, bd, maxiter; warning = false)
Calculate the escape time of a particle `p` in the billiard `bd`, which
is the time until colliding with any "door" in `bd`.
As a "door" is considered any [`FiniteWall`](@ref) with
field `isdoor = true`.

If the particle performs more than `maxiter` collisions without colliding with the
`Door` (i.e. escaping) the returned result is `Inf`.

A warning can be thrown if the result is `Inf`. Enable this using the keyword
`warning = true`.
"""
escapetime(p, bd, t; warning = false) =
    escapetime!(deepcopy(p), bd, t; warning = warning)

function escapetime!(
    p::AbstractParticle{T}, bd::Billiard{T},
    t::Int; warning::Bool=false)::T where {T<:AbstractFloat}

    ei = escapeind(bd)
    if length(ei) == 0
        error("""The billiard does not have any "doors"!""")
    end

    totalt = zero(T); count = zero(t); t_to_write = zero(T)

    while count < t

        i, tmin, pos, vel = bounce!(p, bd)
        t_to_write += tmin

        if typeof(bd[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            totalt += t_to_write
            i ∈ ei &&  break # the collision happens with a Door!

            # set counter
            count += 1
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
    if count ≥ t
        warning && warn("Particle did not escape within max-iter window.")
        return T(Inf)
    end
    return totalt
end



#######################################################################################
## Mean Collision Time
#######################################################################################
function meancollisiontime!(p::AbstractParticle{T}, bd::Billiard{T}, t)::T where {T}

    κ = zero(T)
    ismagnetic = typeof(p) <: MagneticParticle
    ismagnetic && (absω = abs(p.omega))

    count = zero(t); t_to_write = zero(T); colcount = 0

    while count < t

        i, tmin, pos, vel = bounce!(p, bd)

        tmin == Inf && return Inf
        t_to_write += tmin

        if typeof(bd[i]) <: PeriodicWall
            # Pinned particle:
            ismagnetic && t_to_write ≥ 2π/absω && return Inf
            #If not pinned, continue (do not write for PeriodicWall)
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
    meancollisiontime(p, bd, t) -> κ
Compute the mean collision time `κ` of the particle `p` in the billiard `bd` by
evolving for total amount `t` (either float for time or integer for collision number).

Collision times are counted only between obstacles that are *not*
[`PeriodicWall`](@ref).
"""
meancollisiontime(p, bd, t) = meancollisiontime!(deepcopy(p), bd, t)
