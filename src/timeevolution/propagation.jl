export construct, bounce!, evolve, ispinned, evolve!, propagate!

@inline increment_counter(::Int, t_to_write) = 1
@inline increment_counter(::T, t_to_write) where {T<:AbstractFloat} = t_to_write

#####################################################################################
# Bounce
#####################################################################################
"""
    bounce!(p::AbstractParticle, bd::Billiard) → i, t, pos, vel
"Bounce" the particle (advance for one collision) in the billiard.
Takes care of finite-precision issues.

Return:
* index of the obstacle that the particle just collided with
* the time from the previous collision until the current collision `t`
* position and velocity of the particle at the current collision (*after* the
  collision has been resolved!). The position is given in the unit cell of
  periodic billiards. Do `pos += p.current_cell` for the position in real space.

    bounce!(p, bd, raysplit) → i, t, pos, vel
Ray-splitting version of `bounce!`.
"""
@inline function bounce!(p::AbstractParticle{T}, bd::Billiard{T}) where {T}
    i::Int, tmin::T, cp::SV{T} = next_collision(p, bd)
    if tmin != T(Inf)
        o = bd[i]
        relocate!(p, o, tmin, cp)
        resolvecollision!(p, o)
    end
    typeof(p) <: MagneticParticle && (p.center = find_cyclotron(p))
    return i, tmin, p.pos, p.vel
end

#####################################################################################
# Relocate
#####################################################################################
"""
    relocate!(p::AbstractParticle, o::Obstacle, t, cp)
Propagate the particle to `cp` and propagate velocities for time `t`.
Check if it is on the correct side of the obstacle. If not,
change the particle position by [`distance`](@ref) along the [`normalvec`](@ref)
of the obstacle.
"""
@muladd function relocate!(p::AbstractParticle{T},
    o::Obstacle{T}, tmin::T, cp::SV{T}) where {T}

    propagate!(p, cp, tmin) # propagate to collision point
    d = distance(p.pos, o)
    okay = _okay(d, o)
    if !okay
        n = normalvec(o, p.pos)
        p.pos -= d * n
        # due to finite precision this does not always give positive distance
        # but collision takes care of it, as it does not care about collisions
        # very close to the point.
    end
    return okay
end

_okay(d, o::Obstacle) = d ≥ 0.0
_okay(d, o::PeriodicWall) = d ≤ 0.0

"""
    propagate!(p::AbstractParticle, t)
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields.

    propagate!(p, position, t)
Do the same, but take advantage of the already calculated `position` that the
particle should end up at.
"""
@inline function propagate!(p::Particle{T}, t::Real) where {T}
    p.pos += SV{T}(p.vel[1]*t, p.vel[2]*t)
end
@inline propagate!(p::Particle, newpos::SV, t::Real) = (p.pos = newpos)

@inline @muladd function propagate!(p::MagneticParticle{T}, t::Real) where {T}
    ω = p.omega; r = 1/ω
    φ0 = atan(p.vel[2], p.vel[1])
    sinωtφ, cosωtφ = sincos(ω*t + φ0)
    p.pos += SV{T}((sinωtφ - sin(φ0))*r, (-cosωtφ + cos(φ0))*r)
    p.vel = SV{T}(cosωtφ, sinωtφ)
    return
end
@inline @muladd function propagate!(
        p::MagneticParticle{T}, newpos::SV{T}, t) where {T}
    ω = p.omega; φ0 = atan(p.vel[2], p.vel[1])
    p.pos = newpos
    p.vel = SV{T}(cossin(ω*t + φ0))
    return
end


#####################################################################################
# Resolve Collisions
#####################################################################################
"""
    specular!(p::AbstractParticle, o::Obstacle)
Perform specular reflection based on the normal vector of the Obstacle.

In the case where the given obstacle is a `RandomObstacle`, the specular reflection
randomizes the velocity instead (within -π/2+ε to π/2-ε of the normal vector).
"""
@inline function specular!(p::AbstractParticle{T}, o::Obstacle{T})::Nothing where {T}
    n = normalvec(o, p.pos)
    p.vel = p.vel - 2*dot(n, p.vel)*n
    return nothing
end

@inline function specular!(p::AbstractParticle{T}, r::RandomDisk{T})::Nothing where {T}
    n = normalvec(r, p.pos)
    φ = atan(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cos(φ), sin(φ))
    return nothing
end

@inline function specular!(p::AbstractParticle{T}, r::RandomWall{T})::Nothing where {T}
    n = normalvec(r, p.pos)
    φ = atan(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cossin(φ))
    return nothing
end

"""
    periodicity!(p::AbstractParticle, w::PeriodicWall)
Perform periodicity conditions of `w` on `p`.
"""
@inline function periodicity!(p::AbstractParticle, w::PeriodicWall)::Nothing
    p.pos += w.normal
    p.current_cell -= w.normal
    return nothing
end
@inline function periodicity!(p::MagneticParticle, w::PeriodicWall)::Nothing
    p.pos += w.normal
    p.center += w.normal
    p.current_cell -= w.normal
    return nothing
end

"""
    resolvecollision!(p::AbstractParticle, o::Obstacle)
Resolve the collision between particle `p` and obstacle `o`, depending on the
type of `o` (do `specular!` or `periodicity!`).

    resolvecollision!(p, o, T::Function, θ::Function, new_ω::Function)
This is the ray-splitting implementation. The three functions given are drawn from
the ray-splitting dictionary that is passed directly to `evolve!()`. For a calculated
incidence angle φ, if T(φ) > rand(), ray-splitting occurs.
"""
@inline resolvecollision!(p::Particle, o::Obstacle) = specular!(p, o)
@inline resolvecollision!(p::Particle, o::PeriodicWall) = periodicity!(p, o)
@inline resolvecollision!(p::MagneticParticle, o::Obstacle) = specular!(p, o)

resolvecollision!(p::MagneticParticle, o::PeriodicWall) = periodicity!(p, o)

#####################################################################################
# ispinned, Evolve & Construct
#####################################################################################
"""
    ispinned(p::MagneticParticle, bd::Billiard)
Return `true` if the particle is pinned with respect to the billiard.
Pinned particles either have no valid collisions (go in circles forever)
or all their valid collisions are with periodic walls, which again means that
they go in cirles for ever.
"""
ispinned(p::Particle, bd) = false
function _reset_ispinned(p, pos, vel, cc)
    p.pos = pos; p.vel = vel
    p.current_cell = cc; p.center = find_cyclotron(p)
end
function ispinned(p::MagneticParticle, bd::Billiard)
    pos, vel, cc = p.pos, p.vel, p.current_cell
    i, t = bounce!(p, bd)
    if t == Inf; _reset_ispinned(p, pos, vel, cc); return true; end
    if !isperiodic(bd); _reset_ispinned(p, pos, vel, cc); return false; end

    peridx = bd.peridx
    if i ∉ peridx; _reset_ispinned(p, pos, vel, cc); return false; end

    # propagate until 2π/ω
    counter = t; limit = 2π/abs(p.omega)
    while counter ≤ limit
        i, t = bounce!(p, bd)
        if i ∉ peridx; _reset_ispinned(p, pos, vel, cc); return false; end
        counter += t
    end
    _reset_ispinned(p, pos, vel, cc); return true
end

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
        pos = SV(s,-c)/p.ω + prevpos - sc0/p.ω
        vel = SV(s,c)
        push!(poss, pos); push!(vels, vel)
        
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
        pos = prevpos + t*prevvel
        vel = prevvel
        push!(poss, pos); push!(vels, vel)
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

TODO: Update this
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
function construct(p::AbstractParticle{T}, bd::Billiard{T}, t, dt;
                      warning::Bool = true) where {T}

    ts = [zero(T)]
    poss = [p.pos]
    vels = [p.vel]

    count = zero(t)
    t_total = zero(T)
    t_to_write = zero(T)

    prevpos = p.pos + p.current_cell
    prevvel = p.vel

    if ispinned(p, bd)
        warn && @warn "Pinned particle – returning one cycle"
        nposs, nvels, nts = extrapolate(p, prevpos, prevvel, 2π/p.ω, dt)
        append!(ts, nts[2:end] .+ t_total)
        append!(poss, nposs[2:end])
        append!(vels, nvels[2:end])
        return ts, poss, vels
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
                push!(poss, p.pos + p.current_cell)
                push!(vels, p.vel)
            else
                # extrapolate & append
                nposs, nvels, nts = extrapolate(p, prevpos, prevvel, t_to_write, dt)

                append!(ts, nts[2:end] .+ t_total)
                append!(poss, nposs[2:end])
                append!(vels, nvels[2:end]) 
            end

            prevpos = p.pos + p.current_cell
            prevvel = p.vel
            
            t_total += t_to_write
            count += increment_counter(t, t_to_write)
            
            t_to_write = zero(T)
            
        end
        
    end
    
    return ts, poss, vels
end

#reasonable defaults for dt
@inline construct(p::MagneticParticle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    construct(p, bd, t, T(0.01); kwargs...)

@inline construct(p::Particle{T}, bd::Billiard{T}, t; kwargs...) where {T} =
    construct(p, bd, t, T(Inf); kwargs...)

