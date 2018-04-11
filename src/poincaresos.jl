export boundarymap
export psos!, psos

#this function only exists because incidence_angle from raysplitting.jl only works
#if you pass the particle *before* collision, which I cannot do because of bounce!
function reflection_angle(p::AbstractParticle{T}, a::Obstacle{T})::T where {T}
    n = normalvec(a, p.pos)
    inverse_dot = clamp(dot(p.vel, n), -1.0, 1.0)
    φ = acos(inverse_dot)
    if cross2D(p.vel, n) < 0
        φ *= -1
    end
    return φ
end



"""
    arcintervals(bt::Billiard)
Generate an array of `SVector`s, with the `i`th `SVector` containing the arc
length intervals corresponding to the `i`th `Obstacle` in `bt`.

Used by [`boundarymap`](@ref) to compute arc lengths.
"""
function arcintervals(bt::Billiard{T, D}) where {T, D}
    intervals = Vector{SVector{2,T}}(D)
    current = zero(T)
    for i ∈ 1:D
        l = totallength(bt[i]) + current
        intervals[i] = SVector{2,T}(current, l)
        current = l
    end
    return intervals
end



"""
```julia
boundarymap(bt::Billiard, t, ps::Vector{<:AbstractParticle})
boundarymap(bt::Billiard, t, n::Int [, ω])
```
Compute the Poincaré section (also called boundary map) of the
billiard `bt` by evolving each particle for total amount `t` (either float for
time or integer for collision number). See below for the returned values.

If `n::Int` is given instead of `ps`,
then `n` random particles are produced in the billiard table using
[`randominside`](@ref). If `ω` is also given, then the particles are magnetic.

The measurement direction of the arclengths of the individual obstacles
is dictated by `bt`, see [`Billiard`](@ref) for more.

Return
* the arclengths at the collisions `ξs`
* the incidence angles at the collisions `φs`
* obstacle arclength `intervals`

Both `ξs` and `φs` are vectors of `Vector`.
The `i` inner vectors correspond to the results of the `i` initial condition/particle.

The `intervals` is a vector of `SVector`. The `i` entry of `intervals` is the
arclength spanned by the `i` obstacle of the billiard table.
"""
function boundarymap(bt::Billiard{T}, t,
                         ps::Vector{<:AbstractParticle{T}}) where {T}

    params = Vector{T}[]
    angles = Vector{T}[]

    intervals = arcintervals(bt)

    for p ∈ ps
        pparams = T[]
        pangles = T[]
        count = zero(T)
        t_to_write = zero(T)

        while count < t
            i, tmin = bounce!(p,bt)
            t_to_write += tmin

            if typeof(bt[i]) <: PeriodicWall
                continue # do not write output if collision with with PeriodicWall
            else
                push!(pparams, arclength(p, bt[i]) + intervals[i][1])
                push!(pangles, reflection_angle(p, bt[i]))
                # set counter
                count += increment_counter(t, t_to_write)
                t_to_write = zero(T)
            end
        end #time, or collision number, loop

        push!(params, pparams)
        push!(angles, pangles)
    end

    return params, angles, intervals
end

boundarymap(bt::Billiard, t, n::Int) =
    boundarymap(bt, t, [randominside(bt) for i ∈ 1:n])

boundarymap(bt::Billiard, t, n::Int, ω::AbstractFloat) =
    boundarymap(bt, t, [randominside(bt, ω) for i ∈ 1:n])


######################################################################################
# Poincare sos
######################################################################################
psos(p, args...) = psos!(deepcopy(p), args...)

"""
```julia
boundarymap(bt::Billiard, t, ps::Vector{<:AbstractParticle})
boundarymap(bt::Billiard, t, n::Int [, ω])
```

    psos!(p, bt, plane::InfiniteWall, t = 1000) -> poss, vels
Compute the Poincaré section of `p` moving in `bt` when crossing the given `plane`.
Return the positions `poss` and velocities `vels` of the
instances that the crosses the `plane`.

The `plane` can be an [`InfiniteWall`](@ref) of *any* orientation.

Total time of evolution is `t`, which can be either integer or float
(see [`evolve!`](@ref)). Use `psos` if you don't want to mutate the particle.

*Note* - "Pinned" orbits have only 0 or 1 entries in the returned values, depending
on whether they cross the `plane`.
"""
function psos!(
    p::AbstractParticle{T}, bt::Billiard{T}, plane::InfiniteWall{T}, t = 1000) where {T}

    rpos = SV{T}[]
    rvel = SV{T}[]
    count = zero(t)

    while count < t
        # compute collision times
        tmin::T, i::Int = next_collision(p, bt)
        tplane = collisiontime(p, plane)

        # if tplane is smaller, I intersect the section
        if tplane ≥ 0 && tplane < tmin && tplane != Inf
            psos_pos = propagate_pos(p.pos, p, tplane)
            psos_vel = propagate_vel(p, tplane)
            push!(rpos, psos_pos); push!(rvel, psos_vel)
        end

        tmin == Inf && break
        # Now "bounce" the particle normally:
        tmin = relocate!(p, bt[i], tmin)
        resolvecollision!(p, bt[i])

        count += increment_counter(t, tmin)
    end
    return rpos, rvel
end

propagate_vel(p::Particle, t) = p.vel
function propagate_vel(p::MagneticParticle{T}, t) where {T}
    ω = p.omega
    φ0 = atan2(p.vel[2], p.vel[1])
    psos_vel = SV{T}(cos(ω*t + φ0), sin(ω*t + φ0))
end
