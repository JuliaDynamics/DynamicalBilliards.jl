export evolve!, evolve, timeseries!, timeseries

"""
    evolve!([p::AbstractParticle,] bd::Billiard, t)
Evolve the given particle `p` inside the billiard `bd`. If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Return the states of the particle between collisions.

This function mutates the particle, use `evolve` otherwise. If a particle is
not given, a random one is picked through [`randominside`](@ref).

### Returns

* `ct::Vector{T}` : Collision times.
* `poss::Vector{SVector{2,T}}` : Positions at the collisions.
* `vels::Vector{SVector{2,T}})` : Velocities exactly after the collisions.
* `ω`, either `T` or `Vector{T}` : Angular velocity/ies (returned only for magnetic
  particles).

Use the function [`construct`](@ref) to create timeseries of positions and
velocities out of these outputs.

The time `ct[i+1]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting
from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since
`poss[1], vels[1]` are the initial conditions. The angular velocity `ω[i]` is the one
the particle has while propagating from state `poss[i], vels[i]` to `i+1`.

Notice that at any point, the velocity vector `vels[i]` is the one obdained *after*
the specular reflection of the `i-1`th collision.
The function [`construct`](@ref) takes that into account.

### Ray-splitting billiards
    evolve!(p, bd, t, raysplitters)

To implement ray-splitting, the `evolve!` function is supplemented with a
fourth argument, `raysplitters` which is a tuple of [`RaySplitter`](@ref) instances.
For more information and instructions on using ray-splitting
please visit the official documentation.
"""
function evolve!(p::AbstractParticle{T}, bd::Billiard{T}, t;
    warning = false) where {T<:AbstractFloat}

    if t ≤ 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    ismagnetic = typeof(p) <: MagneticParticle
    rt = T[]; push!(rt, 0)
    rpos = SVector{2,T}[]; push!(rpos, p.pos)
    rvel = SVector{2,T}[]; push!(rvel, p.vel)
    ismagnetic && (omegas = T[]; push!(omegas, p.omega); absω = abs(p.omega))

    count = zero(t)
    t_to_write = zero(T)

    if ispinned(p, bd)
        push!(rpos, p.pos)
        push!(rvel, p.vel)
        push!(rt, Inf)
        return (rt, rpos, rvel, p.omega)
    end

    while count < t

        i, tmin, pos, vel = bounce!(p, bd)
        t_to_write += tmin

        if isperiodic(bd) && i ∈ bd.peridx
            continue
        else
            push!(rpos, pos + p.current_cell)
            push!(rvel, vel)
            push!(rt, t_to_write)
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop

    return ismagnetic ? (rt, rpos, rvel, p.omega) : (rt, rpos, rvel)
end

"""
    evolve(p, args...)
Same as [`evolve!`](@ref) but copies the particle instead.
"""
evolve(p::AbstractParticle, args...) = evolve!(copy(p), args...)
evolve(bd::Billiard, args...; kwargs...) =
    evolve!(randominside(bd), bd, args...; kwargs...)



#####################################################################################
# Timeseries
#####################################################################################
"""
    timeseries!(p::AbstractParticle, bd::Billiard, t, [dt])
Evolves the given particle `p` inside the billiard `bd`.  If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Returns the time series for position and velocity as well as the time vector.

The optional argument `dt` is the time step used for interpolating the time
series in between collisions. `dt` is capped by the collision time, as
the interpolation _always_ stops at collisions.
For straight propagation `dt = Inf`, while for magnetic `dt = 0.01`.

For pinned magnetic particles, `timeseries!` issues a warning and returns the
trajectory for one period.

Return:
* x position time-series
* y position time-series
* x velocity time-series
* y velocity time-series
* time vector
"""
function timeseries!(p::AbstractParticle{T}, bd::Billiard{T}, t, dt;
                      warning::Bool = true) where {T}

    ts = [zero(T)]
    x  = [p.pos[1]]; y  = [p.pos[2]]
    vx = [p.vel[1]]; vy = [p.vel[2]]

    count = zero(t)
    t_total = zero(T)
    t_to_write = zero(T)

    prevpos = p.pos + p.current_cell
    prevvel = p.vel

    if ispinned(p, bd)
        warning && @warn "Pinned particle – returning one cycle"
        nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, 2π/p.ω, dt)

        append!(ts, nts[2:end] .+ t_total)
        append!(x, nx[2:end])
        append!(y, ny[2:end])
        append!(vx, nvx[2:end])
        append!(vy, nvy[2:end])

        return x, y, vx, vy, ts
    end

    @inbounds while count < t

        i, ct = bounce!(p, bd)
        t_to_write += ct


        if isperiodic(bd) && i ∈ bd.peridx
            # do nothing at periodic obstacles
            continue
        else
            if t_to_write ≤ dt
                # push collision point only
                push!(ts, t_to_write)
                push!(x, p.pos[1] + p.current_cell[1])
                push!(y, p.pos[2] + p.current_cell[2])

                push!(vx, p.vel[1])
                push!(vy, p.vel[2])
            else
                # extrapolate & append
                nx, ny, nvx, nvy, nts = extrapolate(p, prevpos, prevvel, t_to_write, dt)

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
                     dt::T) where {T <: AbstractFloat}

    φ0 = atan(prevvel[2], prevvel[1])
    s0, c0 = sincos(φ0)

    tvec = 0:dt:ct
    x = Vector{T}(undef, length(tvec))
    y = Vector{T}(undef, length(tvec))
    vx = Vector{T}(undef, length(tvec))
    vy = Vector{T}(undef, length(tvec))

    @inbounds for (i,t) ∈ enumerate(tvec)
        s,c = sincos(p.ω*t + φ0)

        x[i] = s/p.ω + prevpos[1] - s0/p.ω
        y[i] = -c/p.ω + prevpos[2] + c0/p.ω
        vx[i] = s; vy[i] = c
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
                     dt::T) where {T <: AbstractFloat}

    tvec = 0:dt:ct
    poss = Vector{SV{T}}(length(tvec), undef)
    vels = Vector{SV{T}}(length(tvec), undef)

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

#reasonable defaults for dt
@inline timeseries!(p::MagneticParticle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    timeseries!(p, bd, t, T(0.01); kwargs...)

@inline timeseries!(p::Particle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    timeseries!(p, bd, t, T(Inf); kwargs...)


# non-mutating version
"""
    timeseries(p, args...; kwargs...)
Non-mutating version of [`timeseries!`](@ref)
"""
@inline timeseries(p::AbstractParticle, args...; kwargs...) =
    timeseries!(copy(p), args...; kwargs...)
