# ParticlesObstacles.jl must be loaded BEFORE this
export resolvecollision!, collisiontime, propagate!, evolve!, construct, specular!,
periodicity!, propagate_pos, next_collision, escapetime, relocate!

####################################################
## Mathetical/Convenience Functions
####################################################
const sixsqrt = 6sqrt(2)

# Used in relocate:
@inline timeprec(::Type{T}) where {T} = eps(T)^(4/5)

# This timeprec cannot be used for PeriodWall and RaySplitting obstacles with
# MagneticParticle
# because when mangetic and relocating forward you get extremely shallow angles
# and you need huge changes in time for even tiny changes in position
@inline timeprec_forward(::Type{T}) where {T} = eps(T)^(3/4)
@inline timeprec_forward(::Type{BigFloat}) = BigFloat(1e-12)

# Used in check of skip intersection, in `realangle` and collision with Semicircle:
@inline distancecheck(::Type{T}) where {T} = sqrt(eps(T))
@inline distancecheck(::Type{BigFloat}) = BigFloat(1e-8)

"""
    acos1mx(x)
Approximate arccos(1 - x) for x very close to 0.
"""
@inline (acos1mx(x::T)::T) where {T} = sqrt(2x) + sqrt(x)^3/sixsqrt

@inline cross2D(a, b) = a[1]*b[2]-a[2]*b[1]

@inline increment_counter(::Int, t_to_write) = 1
@inline increment_counter(::T, t_to_write) where {T<:AbstractFloat} = t_to_write

####################################################
## Resolve Collisions
####################################################
"""
    specular!(p::AbstractParticle, o::Obstacle)
Perform specular reflection based on the normal vector of the Obstacle.

In the case where the given obstacle is a `RandomObstacle`, the specular reflection
randomizes the velocity instead (within -π/2+ε to π/2-ε of the normal vector).
"""
function specular!(p::AbstractParticle{T}, o::Obstacle{T})::Void where {T}
    n = normalvec(o, p.pos)
    p.vel = p.vel - 2*dot(n, p.vel)*n
    return nothing
end

function specular!(p::AbstractParticle{T}, r::RandomDisk{T})::Void where {T}
    n = normalvec(r, p.pos)
    φ = atan2(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cos(φ), sin(φ))
    return nothing
end

function specular!(p::AbstractParticle{T}, r::RandomWall{T})::Void where {T}
    n = normalvec(r, p.pos)
    φ = atan2(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector{2,T}(cos(φ), sin(φ))
    return nothing
end

"""
    periodicity!(p::AbstractParticle, w::PeriodicWall)
Perform periodicity conditions of `w` on `p`.
"""
function periodicity!(p::AbstractParticle, w::PeriodicWall)::Void
    p.pos += w.normal
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
resolvecollision!(p::AbstractParticle, o::Obstacle)::Void = specular!(p, o)
resolvecollision!(p::AbstractParticle, o::PeriodicWall)::Void =  periodicity!(p, o)

"""
    relocate!(p::AbstractParticle, o::Obstacle, t) -> newt
Propagate the particle's position for time `t`, and check if it is on
the correct side of the obstacle. If not, adjust the time `t` by `timeprec`
and re-evalute until correct. When correct, propagate the particle itself
to the correct position and return the final adjusted time.

Notice that the adjustment is increased geometrically; if one adjustment is not
enough, the adjusted time is multiplied by a factor of 10. This happens as many
times as necessary.
"""
function relocate!(p::AbstractParticle{T}, o::Obstacle{T}, tmin) where {T}
    newpos = propagate_pos(p.pos, p, tmin)
    i = 1
    while distance(newpos, o) < 0
        tmin -= i*timeprec(T)
        newpos = propagate_pos(p.pos, p, tmin)
        i *= 10
    end
    propagate!(p, newpos, tmin)
    return tmin
end

function relocate!(p::AbstractParticle{T}, o::PeriodicWall{T}, tmin) where {T}
    newpos = propagate_pos(p.pos, p, tmin)
    i = 1
    while distance(newpos, o) > 0
        tmin += i*timeprec(T)
        newpos = propagate_pos(p.pos, p, tmin)
        i *= 10
    end
    propagate!(p, newpos, tmin)
    return tmin
end

function relocate!(p::MagneticParticle{T}, o::Obstacle{T}, tmin) where {T}
    newpos = propagate_pos(p.pos, p, tmin)
    i = 1
    while distance(newpos, o) < 0
        tmin -= i*timeprec(T)
        newpos = propagate_pos(p.pos, p, tmin)
        i *= 10
    end
    propagate!(p, newpos, tmin)
    return tmin
end

function relocate!(p::MagneticParticle{T}, o::PeriodicWall{T}, tmin) where {T}
    newpos = propagate_pos(p.pos, p, tmin)
    i = 1
    while distance(newpos, o) > 0
        tmin += i*timeprec_forward(T)
        newpos = propagate_pos(p.pos, p, tmin)
        i *= 10
    end
    propagate!(p, newpos, tmin)
    return tmin
end


# function relocate!(
#     p::MagneticParticle{BigFloat}, o::PeriodicWall{BigFloat}, tmin)::BigFloat
#     newpos = propagate_pos(p.pos, p, tmin)
#     while distance(newpos, o) > 0
#         tmin += timeprec_bigfloatperiodic(T)
#         newpos = propagate_pos(p.pos, p, tmin)
#     end
#     propagate!(p, newpos, tmin)
#     return tmin
# end

####################################################
## Straight Propagation
####################################################
"""
    propagate!(p::AbstractParticle, t)
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields.

For a `Particle` the propagation is a straight line
(i.e. velocity vector is constant). For a `MagneticParticle` the propagation
is circular motion with cyclic frequency `p.omega` and radius `1/p.omega`.

    propagate!(p, position, t)
Do the same, but take advantage of the already calculated `position` that the
particle should end up at.
"""
function propagate!(p::Particle{T}, t::Real) where {T}
    # Set initial conditions
    vx0=p.vel[1]
    vy0=p.vel[2]
    # Set current (final) values for `pos` (`vel` does not change)
    p.pos += SVector{2,T}(vx0*t, vy0*t)
    return
end

function propagate!(p::Particle{T}, newpos::SVector{2,T}, t::Real) where {T}
    p.pos = newpos
    return
end

"""
    propagate_pos(pos, p::Particle{T}, t::Real) where {T}
Perform a "fake" propagation, i.e. propagate a position as if it was the particle's
position.
"""
function propagate_pos(pos, p::Particle{T}, t::Real) where {T}
    vx0=p.vel[1]
    vy0=p.vel[2]
    return SVector{2,T}(pos[1] + vx0*t, pos[2] + vy0*t)
end

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



"""
    next_collision(p, bt) -> (tmin, index)
Return the minimum collision time out of all `collisiontime(p, obst)` for `obst ∈ bt`,
as well as the `index` of the corresponding obstacle.
"""
function next_collision(
    p::AbstractParticle{T}, bt::Vector{<:Obstacle{T}})::Tuple{T,Int} where {T}
    tmin::T = T(Inf)
    ind::Int = 0
    for i in eachindex(bt)
        tcol::T = collisiontime(p, bt[i])
        # Set minimum time:
        if tcol < tmin
            tmin = tcol
            ind = i
        end
    end#obstacle loop
    return tmin, ind
end

function next_collision(
    p::AbstractParticle{T}, bt::Tuple)::Tuple{T,Int} where {T}
    findmin(map(x -> collisiontime(p, x), bt))
end

# testing unrollled:
# using Unrolled
#
# @unroll function min_colt_unrolled(p::AbstractParticle, bt::Tuple)
#     T = eltype(p)
#     tmin = T(Inf)
#     ind = i = 0
#     @unroll for obst in bt
#         i+=1
#         tcol::T = collisiontime(p, obst)
#         # Set minimum time:
#         if tcol < tmin
#             tmin = tcol
#             ind = i
#         end
#     end
#     return tmin, ind
# end

# gives same result as next_collision when benchmarked...

"""
    evolve!(p::AbstractParticle, bt, t [, ray_splitter])
Evolve the given particle `p` inside the billiard table `bt`. If `t` is of type
`AbstractFloat`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.
Return the states of the particle between collisions.

The evolution takes into account the particle's Type.
E.g. if `typeof(p) <: MagneticParticle` then magnetic evolution will take place.

### Return
As noted by the "!" at the end of the function, the call changes
`p` (particle). Most importantly however, this function also returns the main output
expected by a billiard system. This output is a tuple of 3 (or 4) vectors:

* `ct::Vector{T}` : Collision times.
* `poss::Vector{SVector{2,T}}` : Positions during collisions.
* `vels::Vector{SVector{2,T}})` : Velocities **exactly after** the collisions.
* `ω`, either `T` or `Vector{T}` : Angular velocity/ies.

In the case of straight propagation, only the first 3 are returned.
In the case of magnetic propagation, the 4th value is returned as well.
This is either the angular velocity of the particle, or in the case of
ray-splitting it is a vector of the angular velocities at each time step.

The time `ct[i]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting
from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since
`poss[1], vels[1]` are the initial conditions. The angular velocity `ω[i]` is the one
the particle has while propagating from state `poss[i], vels[i]` to `i+1`.

Notice that at any point, the velocity vector `vels[i]` is the one obtained **after**
the specular reflection of the (i-1)th collision.
The function `construct` takes that into account.

### Ray-splitting billiards
To implement ray-splitting, the `evolve!()` function is supplemented with a
fourth argument, `ray_splitter::Dict{Int, Any}`, which maps integers
to some kind of Function container (Tuple or Vector). The functions in this
container are: (φ is the angle of incidence)
* T(φ, pflag, ω) : Transmission probability.
* θ(φ, pflag, ω) : Transmission (aka refraction) angle.
* ω_new(ω, pflag) : Angular velocity after transmission.

For more information and instructions on defining these functions
please visit the official documentation.
"""
function evolve!(p::Particle{T}, bt::Vector{<:Obstacle{T}}, t) where {T<:AbstractFloat}

    if t <= 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    rt = T[]
    rpos = SVector{2,T}[]
    rvel = SVector{2,T}[]
    push!(rpos, p.pos)
    push!(rvel, p.vel)
    push!(rt, zero(T))

    count = zero(T)
    t_to_write = zero(T)

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)

        tmin = relocate!(p, bt[i], tmin)
        t_to_write += tmin

        resolvecollision!(p, bt[i])

        if typeof(bt[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel)
            push!(rt, t_to_write)
            # set counter
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
  return (rt, rpos, rvel)
end


"""
    construct(ct, poss, vels[, ω][, dt=0.01])
Given the main output of this package (see `evolve!()` function) construct the
timeseries of the position and velocity, as well as the time vector.

In case of not given ω (or ω == 0), straight construction takes place.
In case of ω != 0 or ω::Vector magnetic construction takes place.

The additional optional argument of `dt` (only valid for Magnetic construction)
is the timestep with which the timeseries are constructed.

### Return
A tuple of the following:
* xt::Vector{T} : x position time-series
* yt::Vector{T} : y position time-series
* vxt::Vector{T} : x velocity time-series
* vyt::Vector{T} : y velocity time-series
* ts::Vector{T} : time vector
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

####################################################
## Magnetic Propagation
####################################################

function propagate!(p::MagneticParticle{T}, t)::Void where {T}
    # "Initial" conditions
    ω = p.omega
    φ0 = atan2(p.vel[2], p.vel[1])
    # Propagate:
    p.pos += SVector{2, T}(sin(ω*t + φ0)/ω - sin(φ0)/ω,
    -cos(ω*t + φ0)/ω + cos(φ0)/ω)
    p.vel = SVector{2, T}(cos(ω*t + φ0), sin(ω*t + φ0))
    return
end

function propagate!(p::MagneticParticle{T}, newpos::SVector{2,T}, t)::Void where {T}
    # "Initial" conditions
    ω = p.omega
    φ0 = atan2(p.vel[2], p.vel[1])
    # Propagate:
    p.pos = newpos
    p.vel = SVector{2, T}(cos(ω*t + φ0), sin(ω*t + φ0))
    return
end

function propagate_pos(pos, p::MagneticParticle{T}, t) where {T}
    # "Initial" conditions
    ω = p.omega
    vx0= p.vel[1]
    vy0= p.vel[2]
    φ0 = atan2(vy0, vx0)
    # Propagate:
    ppos = SVector{2, T}(sin(ω*t + φ0)/ω - sin(φ0)/ω,
    -cos(ω*t + φ0)/ω + cos(φ0)/ω)
    return pos + ppos
end

"""
    realangle(p::MagneticParticle, o::Obstacle, inter::Vector{SVector}, pc, pr)
Given the intersections `inter` of the trajectory of a magnetic particle `p` with
some obstacle `o`, find which of the two is the "real" one, i.e. occurs first.
Returns the angle of first collision, which is equal to the time to first
collision divided by ω.

The function also takes care of problems that may arise when particles are very
close to the obstacle's boundaries, due to floating-point precision.

(the cyclotron center `pc` and radius `pr` are suplimented for efficiency, since they
have been calculated already)
"""
function realangle(p::MagneticParticle{T}, o::Obstacle{T},
    intersections::Vector{SVector{2, T}}, pc::SVector{2, T}, pr::T)::T where {T}

    ω = p.omega
    P0 = p.pos
    PC = pc - P0
    θ::T = Inf
    for i in intersections
        d2 = dot(i-P0,i-P0) #distance of particle from intersection point
        # Check dot product for close points:
        if d2 ≤ distancecheck(T)
            dotp = dot(p.vel, normalvec(o,  p.pos))
            # Case where velocity points away from obstacle:
            dotp ≥ 0 && continue
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
        # Set minimum angle (first collision)
        if θprime < θ
            θ = θprime
        end
    end
    return θ
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



function evolve!(p::MagneticParticle{T}, bt::Vector{<:Obstacle{T}},
    t; warning::Bool = false) where {T<:AbstractFloat}

    if t <= 0
        throw(ArgumentError("`evolve!()` cannot evolve backwards in time."))
    end

    # if isperiodic(bt) && T == BigFloat
    #     error("Currently periodic+magnetic+BigFloat propagation is not supported :(")
    # end

    ω = p.omega
    absω = abs(ω)
    rt = T[]
    rpos = SVector{2,T}[]
    rvel = SVector{2,T}[]
    push!(rpos, p.pos)
    push!(rvel, p.vel)
    push!(rt, zero(T))

    count = zero(t)
    t_to_write = zero(T)

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)

        if tmin == Inf
            warning && warn("Pinned particle in evolve! (Inf col t)")
            push!(rpos, rpos[end])
            push!(rvel, rvel[end])
            push!(rt, Inf)
            return (rt, rpos, rvel, ω)
        end

        tmin = relocate!(p, bt[i], tmin)
        t_to_write += tmin

        resolvecollision!(p, bt[i])

        # Write output only if the collision was not made with PeriodicWall
        if typeof(bt[i]) <: PeriodicWall
            # Pinned particle:
            if t_to_write ≥ 2π/absω
                warning && warn("Pinned particle in evolve! (completed circle)")
                push!(rpos, rpos[end])
                push!(rvel, rvel[end])
                push!(rt, Inf)
                return (rt, rpos, rvel, ω)
            end
            #If not pinned, continue (do not write for PeriodicWall)
            continue
        else
            push!(rpos, p.pos + p.current_cell)
            push!(rvel, p.vel)
            push!(rt, t_to_write)
            # set counting variable
            count += increment_counter(t, t_to_write)
            t_to_write = zero(T)
        end
    end#time loop
    return (rt, rpos, rvel, ω)
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
        φ0 = atan2(vels[i-1][2], vels[i-1][1])
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
            push!(vxt, cos(ω*td + φ0))
            push!(vyt, sin(ω*td + φ0))
            push!(xt, sin(ω*td + φ0)/ω + x0 - sin(φ0)/ω)  #vy0 is sin(φ0)
            push!(yt, -cos(ω*td + φ0)/ω + y0 + cos(φ0)/ω) #vx0 is cos(φ0)
            push!(ts, t0 + td)
        end#collision time
    end#total time
    return xt, yt, vxt, vyt, ts
end



####################################################
## Escape Times
####################################################
function escapeind(bt)
    j = Int[]
    for (i, obst) in enumerate(bt)
        if typeof(obst) <: FiniteWall && obst.isdoor == true
            push!(j, i)
        end
    end
    return j
end


"""
    escapetime(p, bt, maxiter = 1000; warning = false)
Calculate the escape time of a particle `p` in the billiard table `bt`, which
is the time until colliding with any `Door` in `bt`.
As `Door` is considered any [`FiniteWall`](@ref) with
field `isdoor=true`.

If the particle performs more than `maxiter` collisions without colliding with the
`Door` (i.e. escaping) the returned result is `Inf`.

A warning can be thrown if the result is `Inf`. Enable this using the keyword
`warning = true`.
"""
function escapetime(
    p::AbstractParticle{T}, bt::Vector{<:Obstacle{T}},
    t::Int = 1000; warning::Bool=false)::T where {T<:AbstractFloat}

    ipos = copy(p.pos); ivel = copy(p.vel)
    ei = escapeind(bt)
    if length(ei) == 0
        error("Billiard table does not have any Doors!")
    end

    # Uncomment these to have timeseries outputs
    # rt = T[]
    # rpos = SVector{2,T}[]
    # rvel = SVector{2,T}[]
    # push!(rpos, p.pos)
    # push!(rvel, p.vel)
    # push!(rt, zero(T))

    totalt = zero(T)
    count = zero(t)
    t_to_write = zero(T)

    while count < t
        # Declare these because `bt` is of un-stable type!
        tmin::T, i::Int = next_collision(p, bt)

        tmin = relocate!(p, bt[i], tmin)
        t_to_write += tmin

        resolvecollision!(p, bt[i])

        if typeof(bt[i]) <: PeriodicWall
            continue # do not write output if collision with with PeriodicWall
        else
            # Uncomment these to have timeseries outputs
            # push!(rpos, p.pos + p.current_cell)
            # push!(rvel, p.vel)
            # push!(rt, t_to_write)

            totalt += t_to_write
            i ∈ ei &&  break # the collision happens with a Door!

            # set counter
            count += 1
            t_to_write = zero(T)
        end
    end#time, or collision number, loop
    p.pos = ipos; p.vel = ivel
    if count ≥ t
        warning && warn("Particle did not escape within max-iter window.")
        return Inf
    end
    # return (rt, rpos, rvel)
    return totalt
end
