export evolve!, evolve, construct, extrapolate, timeseries!, timeseries

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


"""
    construct(ct, poss, vels [, ω [, dt=0.01]]) -> xt, yt, vxt, vyt, t
Given the output of [`evolve!`](@ref), construct the
timeseries of the position and velocity, as well as the time vector.
In case of not given `ω` (or `ω == 0`), straight construction takes place.
In case of `ω != 0` or `ω::Vector` magnetic construction takes place.
The additional optional argument of `dt` (only valid for magnetic construction)
is the timestep with which the timeseries are constructed.
Return:
* x position time-series
* y position time-series
* x velocity time-series
* y velocity time-series
* time vector
"""
function construct(t::Vector{T},
    poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}}) where {T}

    xt = [pos[1] for pos in poss]
    yt = [pos[2] for pos in poss]
    vxt = [vel[1] for vel in vels]
    vyt = [vel[2] for vel in vels]
    ct = cumsum(t)
    return xt, yt, vxt, vyt, ct
end

function construct(t::Vector{T},  poss::Vector{SVector{2,T}},
vels::Vector{SVector{2,T}}, ω::T, dt=0.01) where {T}

    dt = T(dt)
    ω == 0 && return construct(t, poss, vels)

    xt = [poss[1][1]]
    yt = [poss[1][2]]
    vxt= [vels[1][1]]
    vyt= [vels[1][2]]
    ts = [t[1]]
    ct = cumsum(t)

    for i in 2:length(t)
        φ0 = atan(vels[i-1][2], vels[i-1][1])
        x0 = poss[i-1][1]; y0 = poss[i-1][2]
        colt=t[i]

        t0 = ct[i-1]
        # Construct proper time-vector
        if colt >= dt
            timevec = collect(0:dt:colt)[2:end]
            timevec[end] == colt || push!(timevec, colt)
        else
            timevec = colt
        end

        for td in timevec
            s, c = sincos(ω*td + φ0)
            s0, c0 = sincos(φ0)
            push!(vxt, c)
            push!(vyt, s)
            push!(xt, s/ω + x0 - s0/ω)  #vy0 is sin(φ0)
            push!(yt, -c/ω + y0 + c0/ω) #vx0 is cos(φ0)
            push!(ts, t0 + td)
        end#collision time
    end#total time
    return xt, yt, vxt, vyt, ts
end



"""
    extrapolate(p, prevpos, prevvel, ct, dt)
one day, this will be a beautiful docstring
"""
function extrapolate(p::MagneticParticle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T) where {T <: AbstractFloat}
    poss = SV{T}[]
    vels = SV{T}[]
    
    φ0 = atan(prevvel[2], prevvel[1])
    s0, c0 = sincos(φ0)
    sc0 = SV(s0, -c0)

    tvec = collect(0:dt:ct)
    
    for t ∈ tvec
        s,c = sincos(p.ω*t + φ0)

        push!(poss, SV(s,-c)/p.ω + prevpos - sc0/p.ω)
        push!(vels, SV(s,c))        
    end

    # finish with ct
    if tvec[end] != ct
        push!(tvec, ct)
        push!(poss, p.pos + p.current_cell)
        push!(vels, p.vel)
    end
    
    return poss, vels, tvec
end

function extrapolate(p::Particle{T}, prevpos::SV{T}, prevvel::SV{T}, ct::T,
                     dt::T) where {T <: AbstractFloat}

    poss = SV{T}[]
    vels = SV{T}[]

    tvec = collect(0:dt:ct)

    for t ∈ tvec
        push!(poss, prevpos + t*prevvel)
        push!(vels, prevvel)
    end
    
    # finish with ct
    if tvec[end] != ct
        push!(tvec, ct)
        push!(poss, p.pos + p.current_cell)
        push!(vels, p.vel)
    end
    
    return poss, vels, tvec
end


"""
    timeseries!(p::AbstractParticle, bd::Billiard, t, [dt])
Evolves the given particle `p` inside the billiard `bd`.  If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Returns the time series for position and velocity as well as the time vector.

The optional argument `dt` is the time step used for interpolating the time 
series in between collisions. For straight propagation, this defaults
 to `Inf`, i.e. only collision points are returned. 
For magnetic propagation, it defaults to 0.01.

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

    ts::Vector{T} = [zero(T)]
    x::Vector{T} = [p.pos[1]]
    y::Vector{T} = [p.pos[2]]
    vx::Vector{T} = [p.vel[1]]
    vy::Vector{T} = [p.vel[2]]

    count = zero(t)
    t_total = zero(T)
    t_to_write = zero(T)

    prevpos = p.pos + p.current_cell
    prevvel = p.vel

    if ispinned(p, bd)
        warn && @warn "Pinned particle – returning one cycle"
        (nx, ny), (nvx, nvy), nts = extrapolate(p, prevpos, prevvel, 2π/p.ω, dt)

        append!(ts, nts[2:end] .+ t_total)
        append!(x, nx[2:end]); append!(y, ny[2:end])
        append!(vx, nvx[2:end]); append!(vy; nvy[2:end])

        return x, y, vx, vy, ts
    end
    
    while count < t

        i, ct = bounce!(p, bd)
        t_to_write += ct

        
        if isperiodic(bd) && i ∈ bd.peridx
            # do nothing at periodic obstacles
            continue
        else
            if t_to_write <= dt
                # push collision point only
                push!(ts, t_to_write)
                push!(x, p.pos[1] + p.current_cell[1])
                push!(y, p.pos[2] + p.current_cell[2])

                push!(vx, p.vel[1])
                push!(vy, p.vel[2])
            else
                # extrapolate & append
                (nx, ny), (nvx, nvy), nts = extrapolate(p, prevpos, prevvel,
                                                        t_to_write, dt)

                append!(ts, nts[2:end] .+ t_total)
                append!(x, nx[2:end]); append!(y, ny[2:end])
                append!(vx, nvx[2:end]); append!(vy; nvy[2:end])                
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

#reasonable defaults for dt
@inline timeseries!(p::MagneticParticle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    construct(p, bd, t, T(0.01); kwargs...)

@inline timeseries!(p::Particle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    construct(p, bd, t, T(Inf); kwargs...)


# non-mutating version
"""
    timeseries(p, args...; kwargs...)
Non-mutating version of [`timeseries!`](@ref)
"""
@inline timeseries(p::AbstractParticle, args...; kwargs...) =
    construct!(copy(p), args...; kwargs...)
