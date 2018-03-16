export collisiontime

"""
    collisiontime(p::AbstractParticle, o::Obstacle)
Calculate the collision time (time-until-collision) between given
particle and obstacle. Returns `Inf` if the collision is not possible *or* if the
collision happens backwards in time.

In the case of magnetic propagation, there are always two possible collisions.
The function internally decides which of the two will occur first, based on the
sign of the angular velocity of the magnetic particle.
"""
function collisiontime(p::Particle{T}, w::Wall{T})::T where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    denom >= 0 ? Inf : dot(w.sp-p.pos, n)/denom
end

function collisiontime(p::Particle{T}, w::FiniteWall{T})::T where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    # case of velocity pointing away of wall:
    denom ≥ 0 && return Inf
    posdot = dot(w.sp-p.pos, n)
    # Case of particle starting behind finite wall:
    posdot ≥ 0 && return Inf
    colt = posdot/denom
    intersection = p.pos .+ colt .* p.vel
    dfc = norm(intersection - w.center)
    if dfc > w.width/2
        return Inf
    else
        return colt
    end
end

function collisiontime(p::Particle{T}, d::Circular{T})::T where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    dotp >=0 && return Inf

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #pointing towards circle center: B < 0
    C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ <= 0 && return Inf
    sqrtD = sqrt(Δ)

    # Closest point:
    t = -B - sqrtD
    t <=0 ? Inf : t
end

function collisiontime(p::Particle{T}, d::Antidot{T})::T where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    if d.pflag == true
        dotp >=0 && return Inf
    end

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B < 0
    C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ <= 0 && return Inf
    sqrtD = sqrt(Δ)

    # Closest point (may be in negative time):
    if dotp < 0
        t = -B - sqrtD
    else
        t = -B + sqrtD
    end

    # If collision time is negative, return Inf:
    t <= 0.0 ? Inf : t
end

function collisiontime(p::Particle{T}, d::Semicircle{T})::T where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ <= 0 && return Inf
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        # these lines make sure that the code works for ANY starting position:
        t = -B - sqrtD
        if t ≤ 0 || distance(p, d) ≤ distancecheck(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel .* t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return Inf
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? Inf : t
end




function collisiontime(p::MagneticParticle{T}, w::Wall{T})::T where {T}
    ω = p.omega
    pc, pr = cyclotron(p)
    P0 = p.pos
    P2P1 = w.ep - w.sp
    P1P3 = w.sp - pc
    # Solve quadratic:
    a = dot(P2P1, P2P1)
    b = 2*dot(P2P1, P1P3)
    c = dot(P1P3, P1P3) - pr^2
    Δ = b^2 -4*a*c
    # Check if line is completely outside (or tangent) of the circle:
    Δ ≤ 0.0 && return Inf
    # Intersection coefficients:
    u1 = (-b - sqrt(Δ))/2a
    u2 = (-b + sqrt(Δ))/2a
    cond1 = (0.0 ≤ u1 ≤ 1.0)
    cond2 = (0.0 ≤ u2 ≤ 1.0)
    # Check if the line is completely inside the circle:
    !cond1 && !cond2 && return Inf
    # Calculate intersection points:
    intersections = SVector{2, T}[]
    cond1 && push!(intersections, w.sp + u1*(w.ep - w.sp))
    cond2 && push!(intersections, w.sp + u2*(w.ep - w.sp))
    # Calculate real time until intersection:
    θ = realangle(p, w, intersections, pc, pr)
    # Collision time, equiv. to arc-length until collision point:
    return θ*pr
end

function collisiontime(p::MagneticParticle{T}, o::Circular{T})::T where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d >= rc + r1) || (d <= abs(rc-r1))
        return Inf
    end
    # Solve quadratic:
    a = (rc^2 - r1^2 + d^2)/2d
    h = sqrt(rc^2 - a^2)
    # Collision points (always 2):
    I1 = SVector{2, T}(
    pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
    pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
    I2 = SVector{2, T}(
    pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
    pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
    # Calculate real time until intersection:
    θ = realangle(p, o, [I1, I2], pc, rc)
    # Collision time, equiv. to arc-length until collision point:
    return θ*rc
end

function collisiontime(p::MagneticParticle{T}, o::Semicircle{T})::T where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d >= rc + r1) || (d <= abs(rc-r1))
        return Inf
    end
    # Solve quadratic:
    a = (rc^2 - r1^2 + d^2)/2d
    h = sqrt(rc^2 - a^2)
    # Collision points (always 2):
    I1 = SVector{2, T}(
    pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
    pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
    I2 = SVector{2, T}(
    pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
    pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
    # Only consider intersections on the "correct" side of Semicircle:
    II = SVector{2,T}[]
    if dot(I1-o.c, o.facedir) < 0 #intersection 1 is OUT
        push!(II, I1)
    end
    if dot(I2-o.c, o.facedir) < 0
        push!(II, I2)
    end
    if length(II) == 0
        return Inf
    end
    # Calculate real time until intersection:
    θ = realangle(p, o, II, pc, rc)
    # Collision time, equiv. to arc-length until collision point:
    return θ*rc
end
