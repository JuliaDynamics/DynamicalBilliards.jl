export bounce!, evolve, ispinned, evolve!, propagate!

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

```julia
bounce!(p, bd, raysplit) → i, t, pos, vel
```
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
