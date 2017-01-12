# ParticlesObstacles.jl must be loaded BEFORE this

####################################################
## Mathetical/Convenience Functions
####################################################
acos1mx(x) = sqrt(2x) + sqrt(x)^3/(6sqrt(2))

####################################################
## Resolve Collisions
####################################################
#All these functions will become dependend on velocity angle for ray-splitting billiards.
"""
```julia
specular!(p::AbstractParticle, o::Obstacle)
```
Perform specular reflection, i.e. set `p.vel = p.vel - 2*dot(n, p.vel)*n` with `n` the
normal vector from the obstacle's boundary.
"""
function specular!{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T})
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
end

"""
    resolvecollision!(p::AbstractParticle, o::Obstacle)
Resolve the collision between particle `p` and obstacle `o`. If the obstacle is not a
periodic wall, the function performs specular reflection. If it is a periodic wall, it
performs the periodicity condition.

`resolvecollision!()` takes special care so that the particle is always inside the correct
side of the billiard table, in order to avoid particle leakage.
Specifically, it calculates the distance from particle and obstacle and,
depending on the obstacle type, makes necessary adjustments by propagating
the particle forwards or backwards in time using **linear** motion.
"""
function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T})

  dist = distance(p, o)
  dt = zero(T)

  if dist < 0.0
    dt = 10lct(p, o, dist) #linearized collision time, should be negative
    # Check the orders of magnitude:
    vxt=p.vel[1]*dt; vyt= p.vel[2]*dt
    #exponent of velocity*time
    ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
    #exponent of position
    pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
    # Ensure that the position will be changed
    if pe - ve > 16.0
      dt *= 10^(pe - ve - 16.0)
    end
    # Propagate backwards:
    p.pos += [p.vel[1]*dt, p.vel[2]*dt]
  end
  # Perform specular reflection:
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
  return dt
end

function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::PeriodicWall{T})

    dist = distance(p, o)
    dt = zero(T)

    if dist > 0.0
      dt = 10lct(p, o, dist) #linearized collision time, should be positive
      # Check the orders of magnitude:
      vxt=p.vel[1]*dt; vyt= p.vel[2]*dt
      #exponent of velocity*timeL
      ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
      #exponent of position
      pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
      # Ensure that the position will be changed
      if pe - ve > 16
        dt *= 10^(pe - ve - 16)
      end
      # Propagate forwards:
      p.pos += [p.vel[1]*dt, p.vel[2]*dt]
    end
  #perform periodicity
  p.pos += o.normal
  p.current_cell -= o.normal
  return dt
end


####################################################
## Linearized collision times
####################################################
"""
```julia
lct(p::AbstractParticle, o::Obstacle, distance)
```
Return the Linearized Collision Time (`-dist/dot(p.vel, normal)`) between particle
and obstacle, given the calculated distance between them.
"""
function lct{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T}, dist::T)
  n = normalvec(o, p.pos)
  t = -dist/dot(p.vel, n)
end

####################################################
## Straight Propagation
####################################################

"""
```julia
collisiontime(p::AbstractParticle, o::Obstacle)
```
Calculate the collision time between given particle and obstacle.

The funtion chooses the appropriate method depending on the type of particle (magnetic or
not) as well as the type of the obstacle. Returns the time that the particle, given its
current position and Type, must be propagated to reach the collision point. This time can
be given directly to `propagate!(p, time)` which brings the particle to the collision point.

In the case of magnetic propagation, there are always two possible collisions. The function
internally decides which of the two will occur first, based on the sign of the angular
velocity of the magnetic particle.
"""
function collisiontime{T<:AbstractFloat}(p::Particle{T}, w::Wall{T})
  n = normalize(w.normal)
  denom = dot(p.vel, n)
  if denom >= 0; return convert(T, Inf); end
  t = dot(w.sp-p.pos, n)/denom
  t <= zero(T) ? convert(T, Inf) : t
end

function collisiontime{T<:AbstractFloat}(p::Particle{T}, d::Disk{T})

  infT = convert(T, Inf)
  dotp = dot(p.vel, normalvec(d, p.pos))
  # Gotta rethink thins for ray spliting
  dotp >=0 && return infT

  dc = p.pos - d.c
  B = dot(p.vel, dc)         #pointing towards circle center: B < 0
  C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
  Δ = B^2 - C

  Δ <= 0 && return infT
  sqrtD = sqrt(Δ)

  # Closest point:
  t = -B - sqrtD
  # Case of being on top of looking inside:
  if t==0.0 && B < 0.0
    t=-2B
  # Case of being inside but closest point is in negative time
  elseif t < 0.0 && C < 0.0
    t = -B + sqrtD
  end

  t <=0 ? convert(T, Inf) : t
end


function collisiontime{T<:AbstractFloat}(p::Particle{T}, d::Circle{T})

  infT = convert(T, Inf)
  dc = p.pos - d.c
  B = dot(p.vel, dc)         #pointing towards circle: B < 0
  C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
  Δ = B^2 - C

  Δ <= 0 && return infT
  sqrtD = sqrt(Δ)

  # Closest point:
  t = -B - sqrtD
  # Case of being on top of looking inside:
  if t==0.0 && B < 0.0
    t=-2B
  # Case of being inside but closest point is in negative time
  elseif t < 0.0 && C < 0.0
    t = -B + sqrtD
  end

  t <=0 ? convert(T, Inf) : t
end


"""
```julia
propagate!(p::AbstractParticle, t)
```
Propagate the particle `p` for given time `t`, changing appropriately the the `p.pos` and
`p.vel` fields.

For a `Particle` the propagation is a straight line (i.e. velocity vector is constant).

For a `MagneticParticle` the propagation is circular motion with cyclic frequency `p.omega`
and radius `1/p.omega`.
"""
function propagate!{T<:AbstractFloat}(p::Particle{T}, t::T)
  # Set initial conditions
  vx0=p.vel[1]
  vy0=p.vel[2]
  # Set current (final) values for `pos` (`vel` is unchanged)
  p.pos += [vx0*t, vy0*t]
end

"""
    evolve!(p::AbstractParticle, bt::Vector{Obstacle}, ttotal)
Evolve the given particle `p` inside the billiard table `bt` for a total amount of time
`ttotal`. Return the states of the particle between collisions.

# Calling
Call the function like:
`t, poss, vels = evolve!(p, bt, ttotal)`
(see "Returns" section for more)

To get the position, velocity and time timeseries from the above output,
use the function `construct`:
`xt, yt, vxt, vyt, ts = construct(t, poss, vels)`
for straight propagation, and:
`xt, yt, vxt, vyt, ts = construct(omega, t, poss, vels)`
for magnetic propagation, with `omega = p.omega` the angular velocity of the particle.

# Returns
    (t::Vector, poss::Vector{SVector{2}}, vels::Vector{SVector{2}})
As noted by the "!" at the end of the function, it changes its argument `p` (particle).
Most importantly however, this function also returns the main output expected by a billiard
system. This output is a tuple of three vectors:
* `t::Vector` : Collision times.
* `poss::Vector{SVector{2}}` : Positions during collisions.
* `vels:: Vector{SVector{2}})` : Velocities **exactly after** the collisions.

The time `t[i]` is the time necessary to reach state `poss[i], vels[i]` starting from the
state `poss[i-1], vels[i-1]`. That is why `t[1]` is always 0 since `poss[0], vels[0]` are
the initial conditions.

Notice that at any point, the velocity vector `vels[i]` is the one obtained **after**
the specular reflection of the (i-1)th collision.
The function `construct` takes that into account.

# Ray-splitting billiards

# Implementation
The function propagates the particle from obstacle to obstacle. At each step the collision
time between all obtacles is calculated and the minimum one is chosen. The particle is
propagated for this amount of time and then the collision between the particle and the
appropriate object is resolved. The loop then begins anew, until the given `ttotal` is
*reached or exceeded* (since only the time between collisions is measured).
"""
function evolve!{T<:AbstractFloat}(p::AbstractParticle{T},
  bt::Vector{Obstacle{T}}, ttotal::T)

  rt = T[]
  rpos = SVector{2,T}[]
  rvel = SVector{2,T}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, zero(T)) #using 0.0 instead works perfectly fine

  infT = convert(T, Inf)
  prev_obst = NullWall(T)
  colobst = NullWall(T)
  tcount = zero(T)
  t_to_write = zero(T)

  while tcount <= ttotal
    tmin = infT

    for obst in bt
      # Always skip collision with previous obstacle if it's a wall
      obst === prev_obst && typeof(prev_obst) <: Wall &&  continue
      tcol = collisiontime(p, obst)
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst = obst
      end
    end#obstacle loop

    propagate!(p, tmin)
    prev_obst = (typeof(colobst) <: PeriodicWall ? colobst.partner : colobst)
    dt = resolvecollision!(p, colobst)
    t_to_write += tmin + dt

    if typeof(colobst) <: PeriodicWall
      continue
    else
      push!(rpos, p.pos + p.current_cell)
      push!(rvel, p.vel)
      push!(rt, t_to_write)
      tcount += t_to_write
      t_to_write = zero(T)
    end
  end#time loop
  return (rt, rpos, rvel)
end

"""
```julia
construct{T<:AbstractFloat}(t::Vector{T},
          poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}})

construct(omega, t, poss, vels, dt=0.01*one(T))
```
Given the main output of this package (see `evolve!()` function) construct the timeseries
of the position and velocity, as well as the time vector.

# Calling
Call this function like:
```julia
xt, yt, vxt, vyt, ts = construct(t, poss, vels)
xt, yt, vxt, vyt, ts = construct(ω, t, poss, vels, dt)
```
for straight and for magnetic propagation respectively.
(`ω = p.omega` the angular velocity of the particle)

The optional argument `dt` settles the frequency points will be written in the output
vectors, which is why its default value depends on whether you have a magnetic or standard
billiard.

# Ray-splitting billiards

# Returns
A tuple of the following:
* xt::Vector{T} : x position time-series
* yt::Vector{T} : y position time-series
* vxt::Vector{T} : x velocity time-series
* vyt::Vector{T} : y velocity time-series
* ts::Vector{T} : Continuous time vector
"""
function construct{T<:AbstractFloat}(t::Vector{T},
  poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}}, dt=0.5*one(T))

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

function propagate!{T<:AbstractFloat}(p::MagneticParticle{T}, t::T)
  # "Initial" conditions
  ω  = p.omega
  vx0= p.vel[1]
  vy0= p.vel[2]
  φ0 = atan2(vy0, vx0)
  # Propagate:
  p.pos += SVector{2, T}(sin(ω*t + φ0)/ω - sin(φ0)/ω, -cos(ω*t + φ0)/ω + cos(φ0)/ω)
  p.vel = SVector{2, T}(cos(ω*t + φ0), sin(ω*t + φ0))
end

"""
```julia
realangle(p::MagneticParticle, o::Obstacle,
          inter::Vector{SVector{2}}, pc::SVector{2}, pr)
```
Given the intersections `inter` of the trajectory of a magnetic particle `p` with some
obstacle `o`, find which of the two is the "real" one, i.e. occurs first.
Returns the angle of first collision.

The function also takes care of problems that may arise when particles are very close to
the obstacle's boundaries, due to Floating-point precision.

(the cyclotron center `pc` and radius `pr` are suplimented for efficiency, since they
have been calculated already)
"""
function realangle{T<:AbstractFloat}(p::MagneticParticle{T}, o::Obstacle{T},
  intersections::Vector{SVector{2, T}}, pc::SVector{2, T}, pr::T)

  ω = p.omega
  P0 = p.pos
  PC = pc - P0
  θ = convert(T, Inf)
  for i in intersections
    d2 = dot(i-P0,i-P0)
    # Check dot product for close points:
    if d2 <= 1e-10
      dotp = dot(p.vel, normalvec(o,  p.pos))
      # Case where velocity points away from obstacle:
      dotp > 0 && continue
    end

    d2r = (d2/(2pr^2))
    # Correct angle value for small arguments:
    θprime = d2r < 1e-16 ? acos1mx(d2r) : acos(1-d2r)

    # Get "side" of i:
    PI = i - P0
    side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
    # Get angle until i (positive number between 0 and π)
    side < 0 && (θprime = abs(2π-θprime))
    # Set minimum angle (first collision)
    if θprime < θ
      θ = θprime
    end
  end
  return θ
end

function collisiontime{T<:AbstractFloat}(p::MagneticParticle{T}, w::Wall{T})
  ω = p.omega
  pc, pr = cyclotron(p)
  P0 = p.pos
  P2P1 = w.ep - w.sp
  P1P3 = w.sp - pc
  # Solve quadratic:
  a::T = dot(P2P1, P2P1)
  b::T = 2*dot(P2P1, P1P3)
  c::T = dot(P1P3, P1P3) - pr^2
  Δ = b^2 -4*a*c
  # Check if line is completely outside (or tangent) of the circle:
  Δ <= zero(T) && return convert(T, Inf)
  # Intersection coefficients:
  u1 = (-b - sqrt(Δ))/2a
  u2 = (-b + sqrt(Δ))/2a
  cond1 = (zero(T) <= u1 <= one(T))
  cond2 = (zero(T) <= u2 <= one(T))
  # Check if the line is completely inside the circle:
  !cond1 && !cond2 && return convert(T, Inf)
  # Calculate intersection points:
  intersections = SVector{2, T}[]
  cond1 && push!(intersections, w.sp + u1*(w.ep - w.sp))
  cond2 && push!(intersections, w.sp + u2*(w.ep - w.sp))

  # Calculate real time until intersection:
  θ = realangle(p, w, intersections, pc, pr)
  # Collision time, equiv. to arc-length until collision point:
  return θ*pr
end

function collisiontime{T<:AbstractFloat}(p::MagneticParticle{T}, o::Circular{T})
  ω = p.omega
  pc, rc = cyclotron(ω, p)
  p1 = o.c
  r1 = o.r
  d = norm(p1-pc)
  if (d >= rc + r1) || (d <= abs(rc-r1))
    return convert(T, Inf)
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


function evolve!{T<:AbstractFloat}(p::MagneticParticle{T},
                bt::Vector{Obstacle{T}}, ttotal::T, dt::T=0.05*one(T))

  ω = p.omega
  rt = T[]
  rpos = SVector{2,T}[]
  rvel = SVector{2,T}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, zero(T))

  infT = convert(T, Inf)
  colobst = NullWall(T)
  tcount = zero(T)
  t_to_write = zero(T)

  while tcount <= ttotal
    tmin = infT

    for obst in bt
      tcol = collisiontime(p, obst)
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst = obst
      end
    end#obstacle loop

    if tmin == infT
      println("pinned particle! (Inf col t)")
      push!(rpos, rpos[end])
      push!(rvel, rvel[end])
      push!(rt, convert(T, Inf))
      return (rt, rpos, rvel)
    end

    propagate!(p, tmin)
    dt = resolvecollision!(p, colobst)
    t_to_write += tmin + dt
    # Write output only if the collision was not made with PeriodicWall
    if typeof(colobst) <: PeriodicWall
      # Pinned particle:
      if t_to_write >= 2π/ω
        println("pinned particle! (completed circle)")
        push!(rpos, rpos[end])
        push!(rvel, rvel[end])
        push!(rt, convert(T, Inf))
        return (rt, rpos, rvel)
      end
      #If not pinned, continue (do not write for PeriodicWall)
      continue
    else
      push!(rpos, p.pos + p.current_cell)
      push!(rvel, p.vel)
      push!(rt, t_to_write)
      tcount += t_to_write
      t_to_write = zero(T)
    end
  end#time loop
  return (rt, rpos, rvel)
end

function construct{T<:AbstractFloat}(ω::T, t::Vector{T},
  poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}}, dt=0.01*one(T))

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
