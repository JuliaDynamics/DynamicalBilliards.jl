export collision, next_collision

#####################################################################################
# Accuracy & Convenience Functions
#####################################################################################
"""
Approximate arccos(1 - x) for x very close to 0.
"""
@inline (acos1mx(x::T)::T) where {T} = sqrt(2x) + sqrt(x)^3/sixsqrt
const sixsqrt = 6sqrt(2)

@inline cross2D(a, b) = a[1]*b[2] - a[2]*b[1]

@inline accuracy(::Type{T}) where {T} = sqrt(eps(T))
@inline accuracy(::Type{BigFloat}) = BigFloat(1e-32)
@inline nocollision(::Type{T}) where {T} = (T(Inf), SV{T}(0.0, 0.0))

#######################################################################################
## Particle
#######################################################################################
"""
    collision(p::AbstractParticle, o::Obstacle) → t, cp
Find the collision (if any) between given particle and obstacle.
Return the time until collision and the estimated collision point `cp`.

Returns `Inf, SV(0, 0)` if the collision is not possible *or* if the
collision happens backwards in time.

**It is the duty of `collision` to avoid incorrect collisions when the particle is
on top of the obstacle (or very close).**
"""
function collision(p::Particle{T}, w::Wall{T}) where {T}
    n = normalvec(w, p.pos)
    denom = dot(p.vel, n)
    if denom ≥ 0.0
        return nocollision(T)
    else
        t = dot(w.sp - p.pos, n)/denom
        return t, p.pos + t * p.vel
    end
end

function collision(p::Particle{T}, w::FiniteWall{T}) where {T}
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

function collision(p::Particle{T}, d::Circular{T}) where {T}

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

function collision(p::Particle{T}, d::Antidot{T}) where {T}

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

function collision(p::Particle{T}, d::Semicircle{T}) where {T}

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


function collision(p::Particle{T}, e::Ellipse{T}) where {T}
    # First check if particle is "looking at" eclipse if it is outside
    if e.pflag
        # These lines may be "not accurate enough" but so far all is good
        dotp = dot(p.vel, normalvec(e, p.pos))
        dotp ≥ 0.0 && return nocollision(T)
    end

    # http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm
    a = e.a; b = e.b
    # Translate particle with ellipse center (so that ellipse lies on [0, 0])
    pc = p.pos - e.c
    # Find μ, ψ for line equation y = μx + ψ describing particle
    μ = p.vel[2]/p.vel[1]
    ψ = pc[2] - μ*pc[1]

    # Determinant and intersection points follow from the link
    denomin = a*a*μ*μ + b*b
    Δ² = denomin - ψ*ψ
    Δ² ≤ 0 && return nocollision(T)
    Δ = sqrt(Δ²); f1 = -a*a*μ*ψ; f2 = b*b*ψ # just factors
    I1 = SV(f1 + a*b*Δ, f2 + a*b*μ*Δ)/denomin
    I2 = SV(f1 - a*b*Δ, f2 - a*b*μ*Δ)/denomin

    d1 = norm(pc - I1); d2 = norm(pc - I2)
    if e.pflag
        return d1 < d2 ? (d1, I1 + e.c) : (d2, I2 + e.c)
    else # inside the ellipse: one collision is _always_ valid
        if d1 < d2
            dmin, Imin = d1, I1
            dmax, Imax = d2, I2
        else
            dmin, Imin = d2, I2
            dmax, Imax = d1, I1
        end

        if dmin < accuracy(T) # special case for being very close to ellipse
            dotp = dot(p.vel, normalvec(e, Imin))
            dotp ≥ 0 && return (dmax, Imax + e.c)
        end
         # check which of the two points is ahead or behind the obstacle
        z1 = dot(pc - Imax, p.vel)
        return z1 < 0 ? (dmax, Imax + e.c) : (dmin, Imin + e.c)
    end
end

#######################################################################################
## Magnetic particle
#######################################################################################
function collision(p::MagneticParticle{T}, w::Wall{T}) where {T}
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

function collision(p::MagneticParticle{T}, o::Circular{T}) where {T}
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
    return θ1 < θ2 ? (θ1*rc, I1) : (θ2*rc, I2)
end

function collision(p::MagneticParticle{T}, o::Semicircle{T}) where {T}
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
    # Only consider intersections on the "correct" side of Semicircle:
    cond1 = dot(I1-o.c, o.facedir) < 0
    cond2 = dot(I2-o.c, o.facedir) < 0
    # Collision time, equiv. to arc-length until collision point:
    θ, I = nocollision(T)
    if cond1 || cond2
        for (Y, cond) in ((I1, cond1), (I2, cond2))
            if cond
                φ = realangle(p, o, Y)
                φ < θ && (θ = φ; I = Y)
            end
        end
    end
    # Collision time = arc-length until collision point
    return θ*rc, I
end

collision(p::MagneticParticle, e::Ellipse) = error(
"Magnetic propagation for Ellipse is not supported :( Consider contributing a "*
"method to `collision(p::MagneticParticle, e::Ellipse)`!")


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
    if d2 ≤ accuracy(T)*accuracy(T)
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
Compute the [`collision`](@ref) across all obstacles in `bd` and find the minimum
one. Return the index of colliding obstacle, the time and the collision point.
"""
function next_collision end

@generated function next_collision(p, bd::Billiard{T, L, BT}) where {T, L, BT}
    out = :(ind = 0; tmin = T(Inf); cp = SV{T}(0.0, 0.0))
    for j=1:L
        push!(out.args, quote
                            let x = bd[$j]
                                tcol, pcol = collision(p, x)
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
