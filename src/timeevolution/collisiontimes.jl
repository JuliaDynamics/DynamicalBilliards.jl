using LinearAlgebra

export collisiontime, realangle
#######################################################################################
## Particle
#######################################################################################
"""
    collisiontime(p::AbstractParticle, o::Obstacle)
Calculate the collision time between given
particle and obstacle. Returns `Inf` if the collision is not possible *or* if the
collision happens backwards in time.

In the case of magnetic propagation, there are always two possible collisions.
The function [`realangle`](@ref) decides which of the two will occur first,
based on the sign of the angular velocity of the magnetic particle.

**Warning** - It is the duty of the `collisiontime` function to *not* lead to
double collisions (e.g. avoid collisions when the particle is already on the obstacle)!
"""
function collisiontime(p::Particle{T}, w::Wall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    denom >= 0.0 ? T(Inf) : dot(w.sp-p.pos, n)/denom
end

function collisiontime(p::Particle{T}, w::FiniteWall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    # case of velocity pointing away of wall:
    denom ≥ 0.0 && return Inf
    posdot = dot(w.sp-p.pos, n)
    # Case of particle starting behind finite wall:
    posdot ≥ 0.0 && return Inf
    colt = posdot/denom
    intersection = p.pos + colt * p.vel
    dfc = norm(intersection - w.center)
    if dfc > w.width/2
        return T(Inf)
    else
        return colt
    end
end

@muladd function collisiontime(p::Particle{T}, d::Circular{T}) where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    dotp >= 0.0 && return T(Inf)

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #pointing towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B*B - C

    Δ <= 0.0 && return T(Inf)
    sqrtD = sqrt(Δ)

    # Closest point:
    t = -B - sqrtD
    t <= 0.0 ? T(Inf) : t
end

@muladd function collisiontime(p::Particle{T}, d::Antidot{T})::T where {T}

    dotp = dot(p.vel, normalvec(d, p.pos))
    if d.pflag == true
        dotp ≥ 0 && return T(Inf)
    end

    dc = p.pos - d.c
    B = dot(p.vel, dc)           #velocity towards circle center: B < 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ <= 0 && return T(Inf)
    sqrtD = sqrt(Δ)

    # Closest point (may be in negative time):
    if dotp < 0
        t = -B - sqrtD
    else
        t = -B + sqrtD
    end

    # If collision time is negative, return Inf:
    t <= 0.0 ? T(Inf) : t
end

@muladd function collisiontime(p::Particle{T}, d::Semicircle{T})::T where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
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
    newpos = p.pos + p.vel * t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return Inf
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? Inf : t
end

@muladd function collisiontime(p::Particle{T}, e::Ellipse{T})::T where {T}
    # First check if particle is "looking at" eclipse if it is outside
    if e.pflag
        # These lines may be "not accurate enough" but so far all is good
        dotp = dot(p.vel, normalvec(e, p.pos))
        dotp >= 0.0 && return T(Inf)
    end

    # http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm
    a = e.a; b = e.b
    # Translate particle with ellipse center (so that ellipse lies on [0, 0])
    pc = p.pos - e.c
    # Find μ, ψ for line equation y = μx + ψ describing particle
    μ = p.vel[2]/p.vel[1]
    ψ = pc[2] - μ*pc[1]

    # Determinant and intersection points follow from the link
    denominator = a*a*μ*μ + b*b
    Δ² = denominator - ψ*ψ
    Δ² ≤ 0 && return T(Inf)
    Δ = sqrt(Δ²); f1 = -a*a*μ*ψ; f2 = b*b*ψ # just factors
    I1 = SV(f1 + a*b*Δ, f2 + a*b*μ*Δ)/denominator
    I2 = SV(f1 - a*b*Δ, f2 - a*b*μ*Δ)/denominator

    d1 = norm(pc - I1); d2 = norm(pc - I2)
    if e.pflag
        # There HAS to be a way to KNOW which of the two I1, I2 is the "real closest"
        # one. Analytically I mean. But I haven't found yet.
        return min(d1, d2)
    else # If inside the ellipse, only one collision valid
        if min(d1, d2) < distancecheck(T)
            return max(d1,d2)
        else
            dotp = dot(p.vel, (I1 - pc))
            return dotp ≥ 0.0 ? d1 : d2
        end
    end
end

#######################################################################################
## Magnetic particle
#######################################################################################
@muladd function collisiontime(p::MagneticParticle{T}, w::Wall{T})::T where {T}
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
    Δ ≤ 0.0 && return Inf
    # Intersection coefficients:
    u1 = (-b - sqrt(Δ))/2a
    u2 = (-b + sqrt(Δ))/2a
    cond1 = 0.0 ≤ u1 ≤ 1.0
    cond2 = 0.0 ≤ u2 ≤ 1.0
    # Check if the line is completely inside the circle:
    if cond1 || cond2
        # Calculate real angle until intersection:
        θ1 = cond1 ? (I1 = w.sp + u1*(w.ep-w.sp); realangle(p, w, I1)) : T(Inf)
        θ2 = cond2 ? (I2 = w.sp + u2*(w.ep-w.sp); realangle(p, w, I2)) : T(Inf)
        # Collision time, equiv. to arc-length until collision point:
        return min(θ1, θ2)*pr
    else
        return Inf
    end
end

@muladd function collisiontime(p::MagneticParticle{T}, o::Circular{T})::T where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d ≥ rc + r1) || (d ≤ abs(rc-r1))
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
    if typeof(o) <: Antidot && o.pflag == false
        d1 = norm(p.pos - I1); d2 = norm(p.pos - I2)
        j = d1 < d2 ? 2 : 1
        if min(d1, d2) ≤ distancecheck(T)
            return realangle(p, o, (I1, I2)[j])
        end
    end
    # Calculate real time until intersection:
    θ1 = realangle(p, o, I1)
    θ2 = realangle(p, o, I2)
    # Collision time, equiv. to arc-length until collision point:
    return min(θ1, θ2)*rc
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
    if d2 ≤ distancecheck(T)
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
    side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
    # Get angle until i (positive number between 0 and 2π)
    side < 0 && (θprime = T(2π-θprime))
    return θprime
end



#######################################################################################
## next_collision
#######################################################################################
"""
    next_collision(p::AbstractParticle, bd::Billiard) -> tmin, i
Compute the [`collisiontime`](@ref) across all obstacles in `bd`, find the minimum
one and return this time as well as the index of the obstacle that the time
corresponds to.
"""
function next_collision end

@generated function next_collision(p, bd::Billiard{T, L, BT}) where {T, L, BT}
    out = :(i = 0; ind = 0; tmin = T(Inf))
    for j=1:L
        push!(out.args, quote
                            let x = bd[$j]
                                tcol = collisiontime(p, x)
                                # Set minimum time:
                                if tcol < tmin
                                  tmin = tcol
                                  ind = $j
                                end
                            end
                        end
                        )
    end
    push!(out.args, :(return tmin, ind))
    return out
end
