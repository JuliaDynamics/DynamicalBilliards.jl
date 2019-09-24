# This file must be loaded _after_ raysplitting

export evolve!, evolve, timeseries!, timeseries
export visited_obstacles, visited_obstacles!

"""
    evolve!([p::AbstractParticle,] bd::Billiard, t)
Evolve the given particle `p` inside the billiard `bd`. If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Return the states of the particle between collisions.

This function mutates the particle, use `evolve` otherwise. If a particle is
not given, a random one is picked through [`randominside`](@ref).

### Return

* `ct::Vector{T}` : Collision times.
* `poss::Vector{SVector{2,T}}` : Positions at the collisions.
* `vels::Vector{SVector{2,T}})` : Velocities exactly after the collisions.
* `ω`, either `T` or `Vector{T}` : Angular velocity/ies (returned only for magnetic
  particles).

The time `ct[i+1]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting
from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since
`poss[1], vels[1]` are the initial conditions. The angular velocity `ω[i]` is the one
the particle has while propagating from state `poss[i], vels[i]` to `i+1`.

Notice that at any point, the velocity vector `vels[i]` is the one obdained *after*
the specular reflection of the `i-1`th collision.

### Ray-splitting billiards
    evolve!(p, bd, t, raysplitters)

To implement ray-splitting, the `evolve!` function is supplemented with a
fourth argument, `raysplitters` which is a tuple of [`RaySplitter`](@ref) instances.
Notice that `evolve` **always mutates the billiard** if ray-splitting is used!
For more information and instructions on using ray-splitting
please visit the official documentation.
"""
function evolve!(p::AbstractParticle{T}, bd::Billiard{T}, t, raysplitters = nothing;
    warning = false) where {T<:AbstractFloat}

    if t ≤ 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end
    if ispinned(p, bd)
        push!(rpos, p.pos)
        push!(rvel, p.vel)
        push!(rt, Inf)
        return (rt, rpos, rvel, p.ω)
    end

    ismagnetic = p isa MagneticParticle
    isray = !isa(raysplitters, Nothing)
    isray && acceptable_raysplitter(raysplitters, bd)
    raysidx = raysplit_indices(bd, raysplitters)
    ismagnetic && isray && (omegas = [p.ω])

    rt = T[0.0]; rpos = [p.pos]; rvel = [p.vel]
    count = zero(t)
    t_to_write = zero(T)
    if typeof(t) == Int
        for zzz in (rt, rpos, rvel)
            sizehint!(zzz, t)
        end
    end

    while count < t

        i, tmin, pos, vel = bounce!(p, bd, raysidx, raysplitters)
        t_to_write += tmin

        if isperiodic(i, bd)
            continue
        else
            push!(rpos, pos + p.current_cell)
            push!(rvel, vel)
            push!(rt, t_to_write)
            ismagnetic && isray && push!(omegas, p.ω)
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop

    # Return stuff
    if ismagnetic && isray
        return (rt, rpos, rvel, omegas)
    elseif ismagnetic
        return (rt, rpos, rvel, p.ω)
    else
        return (rt, rpos, rvel)
    end
end

"""
    evolve(p, args...)
Same as [`evolve!`](@ref) but copies the particle instead.
"""
evolve(p::AbstractParticle, args...) = evolve!(copy(p), args...)
evolve(bd::Billiard, args...; kwargs...) =
    evolve!(randominside(bd), bd, args...; kwargs...)

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

#####################################################################################
# Timeseries
#####################################################################################
"""
    timeseries!([p::AbstractParticle,] bd::Billiard, t; dt, warning)
Evolves the given particle `p` inside the billiard `bd`.  If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Returns the time series for position and velocity as well as the time vector.

This function mutates the particle, use `timeseries` otherwise. If a particle is
not given, a random one is picked through [`randominside`](@ref).

The keyword argument `dt` is the time step used for interpolating the time
series in between collisions. `dt` is capped by the collision time, as
the interpolation _always_ stops at collisions.
For straight propagation `dt = Inf`, while for magnetic `dt = 0.01`.

For pinned magnetic particles, `timeseries!` issues a warning and returns the
trajectory of the particle. If `t` is integer, the trajectory is evolved for
one full circle only

Return:
* x position time-series
* y position time-series
* x velocity time-series
* y velocity time-series
* time vector

### Ray-splitting billiards
    timeseries!(p, bd, t, raysplitters; ...)

To implement ray-splitting, the `timeseries!` function is supplemented with a
fourth argument, `raysplitters` which is a tuple of [`RaySplitter`](@ref) instances.
Notice that `timeseries` **always mutates the billiard** if ray-splitting is used!
For more information and instructions on using ray-splitting
please visit the official documentation.
"""
function timeseries!(p::AbstractParticle{T}, bd::Billiard{T}, t, raysplitters = nothing;
                     dt = typeof(p) <: Particle ? T(Inf) : T(0.01),
                     warning::Bool = true) where {T}

    ts = [zero(T)]
    x  = [p.pos[1]]; y  = [p.pos[2]]
    vx = [p.vel[1]]; vy = [p.vel[2]]

    ismagnetic = p isa MagneticParticle
    prevω = ismagnetic ? p.ω : T(0)
    isray = !isa(raysplitters, Nothing)
    isray && acceptable_raysplitter(raysplitters, bd)
    raysidx = raysplit_indices(bd, raysplitters)

    count = zero(t)
    t_total = zero(T)
    t_to_write = zero(T)

    prevpos = p.pos + p.current_cell
    prevvel = p.vel

    if ispinned(p, bd)
        warning && @warn "Pinned particle detected!"

        #return one cycle if t is a collision number
        t_ret = typeof(t) <: Integer ? 2π/p.ω : T(t)
        nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, t_ret, dt, p.ω)

        append!(ts, nts[2:end] .+ t_total)
        append!(x, nx[2:end])
        append!(y, ny[2:end])
        append!(vx, nvx[2:end])
        append!(vy, nvy[2:end])

        return x, y, vx, vy, ts
    end

    @inbounds while count < t

        ismagnetic && isray && (prevω = p.ω)
        i, ct = bounce!(p, bd, raysidx, raysplitters)
        t_to_write += ct

        if isperiodic(i, bd)
            # do nothing at periodic obstacles
            continue
        else
            if t_to_write ≤ dt
                # push collision point only
                push!(ts, t_to_write + t_total)
                push!(x, p.pos[1] + p.current_cell[1])
                push!(y, p.pos[2] + p.current_cell[2])

                push!(vx, p.vel[1])
                push!(vy, p.vel[2])
            else
                # extrapolate & append
                nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, t_to_write, dt, prevω)

                append!(ts, nts[2:end] .+ t_total)
                append!(x, nx[2:end])
                append!(y, ny[2:end])
                append!(vx, nvx[2:end])
                append!(vy, nvy[2:end])
            end

            prevpos = p.pos + p.current_cell
            prevvel = p.vel

            t_total += t_to_write
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)

        end
    end
    return x, y, vx, vy, ts
end

function extrapolate(p::MagneticParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T, ω::T) where {T <: AbstractFloat}

    φ0 = atan(prevvel[2], prevvel[1])
    s0, c0 = sincos(φ0)

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        s,c = sincos(ω*t + φ0)
        x[i] = s/ω + prevpos[1] - s0/ω
        y[i] = -c/ω + prevpos[2] + c0/ω
        vx[i] = c; vy[i] = s
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        push!(x, p.pos[1] + p.current_cell[1])
        push!(y, p.pos[2] + p.current_cell[2])
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

function extrapolate(p::Particle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T, ω::T) where {T <: AbstractFloat}

    tvec = collect(0:dt:ct)
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        x[i] = prevpos[1] + t*prevvel[1]
        y[i] = prevpos[2] + t*prevvel[2]
        vx[i], vy[i] = prevvel
    end

    # finish with ct
    @inbounds if tvec[end] != ct
        push!(tvec, ct)
        push!(x, p.pos[1] + p.current_cell[1])
        push!(y, p.pos[2] + p.current_cell[2])
        push!(vx, p.vel[1]); push!(vy, p.vel[2])
    end

    return x, y, vx, vy, tvec
end

# non-mutating version
"""
    timeseries(p, args...; kwargs...)
Non-mutating version of [`timeseries!`](@ref)
"""
@inline timeseries(p::AbstractParticle, args...; kwargs...) =
    timeseries!(copy(p), args...; kwargs...)
