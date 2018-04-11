export psoscut!, psoscut

psoscut(p, args...) = psoscut!(deepcopy(p), args...)

"""
    psoscut!(p, bt, plane::InfiniteWall, t = 1000) -> poss, vels
Compute the Poincaré section of `p` moving in `bt` when crossing the given `plane`.
Return the positions `poss` and velocities `vels` of the
instances that the crosses the `plane`.

The `plane` can be an [`InfiniteWall`](@ref) of *any* orientation.

Total time of evolution is `t`, which can be either integer or float
(see [`evolve!`](@ref)). Use `psoscut` if you don't want to mutate the particle.
"""
function psoscut!(
    p::AbstractParticle{T}, bt::Billiard{T}, plane::InfiniteWall{T}, t = 1000) where {T}

    rpos = SV{T}[]
    rvel = SV{T}[]
    count = zero(t)

    while count < t
        # compute collision times
        tmin::T, i::Int = next_collision(p, bt)
        tplane = collisiontime(p, plane)

        # if tplane is smaller, I intersect the section
        if tplane >= 0 && tplane < tmin
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
