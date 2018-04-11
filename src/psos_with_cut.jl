export psoscut!, psoscut

psoscut(p, args...) = psoscut!(deepcopy(p), args...)

function psoscut!(
    p::AbstractParticle{T}, bt::Billiard{T}, t, plane::InfiniteWall{T}) where {T}

    rpos = SVector{2,T}[]
    rvel = SVector{2,T}[]
    count = zero(t)

    while count < t
        # compute collision times
        tmin::T, i::Int = next_collision(p, bt)
        tplane = collisiontime(p, plane)

        # if tplane is smaller, I intersect the section
        if tplane < tmin
            psos_pos = propagate_pos(p.pos, p, tplane)
            psos_vel = propagate_vel(p, tplane)
            push!(rpos, psos_pos); push!(rvel, psos_vel)
        end

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
