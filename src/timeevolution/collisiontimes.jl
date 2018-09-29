export collisiontime

#####################################################################################
# Accuracy & Convenience Functions
#####################################################################################
"""
Approximate arccos(1 - x) for x very close to 0.
"""
@inline (acos1mx(x::T)::T) where {T} = sqrt(2x) + sqrt(x)^3/sixsqrt
const sixsqrt = 6sqrt(2)

@inline cross2D(a, b) = a[1]*b[2] - a[2]*b[1]

@inline accuracy(::Type{T}) where {T} = eps(T)^(3/4)
@inline accuracy(::Type{BigFloat}) = BigFloat(1e-32)
@inline nocollision(::Type{T}) where {T} = (T(Inf), SV{T}(0.0, 0.0))

#######################################################################################
## Particle
#######################################################################################
"""
    collisiontime(p::AbstractParticle, o::Obstacle) -> t, cp
Calculate the collision time between given
particle and obstacle. Return the time and the estimated collision point `cp`.

Returns `Inf, SV(Inf, Inf)` if the collision is not possible *or* if the
collision happens backwards in time.

**It is the duty of `collisiontime` to avoid incorrect collisions when the particle is
on top of the obstacle (or extremely close).**
"""
@muladd function collisiontime(p::Particle{T}, w::Wall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision(T)
    else
        t = dot(w.sp - p.pos, n)/denom
        return t, p.pos + t * p.vel
    end
end

function collisiontime(p::Particle{T}, w::FiniteWall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    # case of velocity pointing away of wall:
    denom ≥ 0.0 && return nocollision(T)
    posdot = dot(w.sp-p.pos, n)
    # Case of particle starting behind finite wall:
    posdot ≥ 0.0 && return nocollision(T)
    colt = posdot/denom
    i = p.pos + colt * p.vel
    dfc = norm(i - w.center)
    if dfc > w.width/2
        return nocollision(T)
    else
        return colt, i
    end
end

@muladd function collisiontime(p::Particle{T}, d::Circular{T}) where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    dotp ≥ 0.0 && return nocollision(T)

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #pointing towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B*B - C

    Δ ≤ 0.0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    # Closest point:
    t = -B - sqrtD
    return t, p.pos + t * p.vel
end

@muladd function collisiontime(p::Particle{T}, d::Antidot{T}) where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    if d.pflag == true
        dotp ≥ 0 && return nocollision(T)
    end

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #velocity towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    # Closest point (may be in negative time):
    if dotp < 0
        t = -B - sqrtD
    else
        t = -B + sqrtD
    end

    # If collision time is negative, return Inf:
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t * p.vel)
end

@muladd function collisiontime(p::Particle{T}, d::Semicircle{T}) where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        t = -B - sqrtD
        # these lines make sure that the code works for ANY starting position:
        if t ≤ 0 || distance(p, d) ≤ accuracy(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel * t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return nocollision(T)
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t*p.vel)
end



#######################################################################################
## Magnetic particle
#######################################################################################
@muladd function collisiontime(p::MagneticParticle{T}, w::Wall{T}) where {T}
    ω = p.omega
    pc, pr = cyclotron(p)
    P0 = p.pos
    P2P1 = w.ep - w.sp
    P1P3 = w.sp - pc
    # Solve quadratic:
    a = dot(P2P1, P2P1)
    b = 2*dot(P2P1, P1P3)
    c = dot(P1P3, P1P3) - pr*pr
    Δ = b^2 -4*a*c
    # Check if line is completely outside (or tangent) of the circle:
    Δ ≤ 0.0 && return nocollision(T)
    # Intersection coefficients:
    u1 = (-b - sqrt(Δ))/2a
    u2 = (-b + sqrt(Δ))/2a
    cond1 = 0.0 ≤ u1 ≤ 1.0
    cond2 = 0.0 ≤ u2 ≤ 1.0
    # Check if the line (wall) is completely inside the circle:
    θ, I = nocollision(T)
    if cond1 || cond2
        dw = w.ep - w.sp
        for (u, cond) in ((u1, cond1), (u2, cond2))
            Y =  w.sp + u*dw
            if cond
                φ = realangle(p, w, Y)
                φ < θ && (θ = φ; I = Y)
            end
        end
    end
    # Collision time = arc-length until collision point
    return θ*pr, I
end

@muladd function collisiontime(p::MagneticParticle{T}, o::Circular{T}) where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d >= rc + r1) || (d <= abs(rc-r1))
        return nocollision(T)
    end
    # Solve quadratic:
    a = (rc^2 - r1^2 + d^2)/2d
    h = sqrt(rc^2 - a^2)
    # Collision points (always 2):
    I1 = SVector{2, T}(
        pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
        pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d
    )
    I2 = SVector{2, T}(
        pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
        pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d
    )
    # Calculate real time until intersection:
    θ1 = realangle(p, o, I1)
    θ2 = realangle(p, o, I2)
    # Collision time, equiv. to arc-length until collision point:
    j = θ1 < θ2 ? 1 : 2
    return (θ1, θ2)[j]*rc, (I1, I2)[j]
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
        pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d
    )
    I2 = SVector{2, T}(
        pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
        pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d
    )
    # Only consider intersections on the "correct" side of Semicircle:
    cond1 = dot(I1-o.c, o.facedir) < 0
    cond2 = dot(I2-o.c, o.facedir) < 0
    if cond1 || cond2
        # Calculate real angle until intersection:
        θ1 = cond1 ? realangle(p, o, I1) : T(Inf)
        θ2 = cond2 ? realangle(p, o, I2) : T(Inf)
        # Collision time, equiv. to arc-length until collision point:
        return min(θ1, θ2)*rc
    else
        return Inf
    end
end

"""
    realangle(p::MagneticParticle, o::Obstacle, I) -> θ
Given the intersection point `I` of the trajectory of a magnetic particle `p` with
some obstacle `o`, find the real angle that will be spanned until the particle
collides with the obstacle.

The function also takes care of problems that may arise when particles are very
close to the obstacle's boundaries, due to floating-point precision.
"""
function realangle(p::MagneticParticle{T}, o::Obstacle{T}, i::SV{T})::T where {T}

    pc = p.center; pr = p.r; ω = p.omega
    P0 = p.pos
    PC = pc - P0
    d2 = dot(i-P0,i-P0) #distance of particle from intersection point
    # Check dot product for close points:
    if d2 ≤ accuracy(T)
        dotp = dot(p.vel, normalvec(o,  p.pos))
        # Case where velocity points away from obstacle:
        dotp ≥ 0 && return T(Inf)
    end

    d2r = (d2/(2pr^2))
    d2r > 2 && (d2r = T(2.0))
    # Correct angle value for small arguments (between 0 and π/2):
    θprime = d2r < 1e-3 ? acos1mx(d2r) : acos(1.0 - d2r)

    # Get "side" of i:
    PI = i - P0
    side = cross2D(PI, PC)*ω
    # Get angle until i (positive number between 0 and 2π)
    side < 0 && (θprime = T(2π-θprime))
    return θprime
end



#######################################################################################
## next_collision
#######################################################################################
"""
    next_collision(p::AbstractParticle, bd::Billiard) -> i, tmin, cp
Compute the [`collisiontime`](@ref) across all obstacles in `bd` and find the minimum
one. Return the index of colliding obstacle, the time and the collision point.
"""
function next_collision end

@generated function next_collision(p, bd::Billiard{T, L, BT}) where {T, L, BT}
    out = :(ind = 0; tmin = T(Inf); cp = SV{T}(0.0, 0.0))
    for j=1:L
        push!(out.args, quote
                            let x = bd[$j]
                                tcol, pcol = collisiontime(p, x)
                                # Set minimum time:
                                if tcol < tmin
                                  tmin = tcol
                                  ind = $j
                                  cp = pcol
                                end
                            end
                        end
                        )
    end
    push!(out.args, :(return ind, tmin, cp))
    return out
end
