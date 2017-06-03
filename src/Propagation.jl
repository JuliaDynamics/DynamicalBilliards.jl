# ParticlesObstacles.jl must be loaded BEFORE this
export resolvecollision!, collisiontime, propagate!, evolve!, construct, specular!,
periodicity!

####################################################
## Mathetical/Convenience Functions
####################################################
"""
    acos1mx(x)
Approximate arccos(1 - x) for x very close to 0.
"""
acos1mx(x) = sqrt(2x) + sqrt(x)^3/(6sqrt(2))

"""
```julia
cross2D(a, b) = a[1]*b[2]-a[2]*b[1]
```
Return the 3rd element of the cross product of `a` and `b`, assuming their
extention from 2D to 3D arrays simply pushes a 0 at their ends.

Equivalent to the determinant of the matrix [a b].
"""
cross2D{T<:Real, P<:Real}(a::SVector{2,T}, b::SVector{2,P}) = a[1]*b[2]-a[2]*b[1]
cross2D(a::AbstractArray, b::AbstractArray) = a[1]*b[2]-a[2]*b[1]

increment_counter(::Int, t_to_write)::Int = 1
increment_counter(::Float64, t_to_write)::Float64 = t_to_write
####################################################
## Resolve Collisions
####################################################
"""
```julia
relocate!(p::AbstractParticle, o::Obstacle, distance)
```
Relocate the particle `p` with respect to the obstacle `o` along the direction of
the normal vector of `o`, so that it is on the correct side of the obstacle.

Internally takes care of problems of finite accuracy of `Float64`.
"""
function relocate!(p::AbstractParticle, o::Obstacle, dist::Float64)::Void
  n = normalvec(o, p.pos)
  ex = Int(floor(exponent(dist)*0.3010299956639812))
  if ex < -15
    p.pos -= 10^(-ex - 15)*2*dist*n
  else
    p.pos -= 2*dist*n
  end
  newdist = distance(p, o)
  if sign(newdist) == sign(dist)
    error("New distance should have opposite sign from old!")
  end
  return
end

"""
```julia
specular!(p::AbstractParticle, o::Obstacle)
```
Perform specular reflection based on the normal vector of the Obstacle.

In the case where the given obstacle is a `RandomObstacle`, the specular reflection
randomizes the velocity instead (within -π/2 to π/2 of the normal vector).
"""
function specular!(p::AbstractParticle, o::Obstacle)
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
end

function specular!(p::AbstractParticle, r::RandomDisk)
    n = normalvec(r, p.pos)
    φ = atan2(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector(cos(φ), sin(φ))
end

function specular!(p::AbstractParticle, r::RandomWall)
    n = normalvec(r, p.pos)
    φ = atan2(n[2], n[1]) + 0.95(π*rand() - π/2) #this cannot be exactly π/2
    p.vel = SVector(cos(φ), sin(φ))
end

"""
```julia
periodicity!(p::AbstractParticle, w::PeriodicWall)
```
Perform periodicity conditions of `w` on `p`.
"""
function periodicity!(p::AbstractParticle, w::PeriodicWall)
  p.pos += w.normal
  p.current_cell -= w.normal
end

"""
    resolvecollision!(p::AbstractParticle, o::Obstacle)
Resolve the collision between particle `p` and obstacle `o`, depending on the
type of `o`. `resolvecollision!()` takes special care so that the particle is always
inside the correct side of the billiard table, in order to avoid particle leakage.

    resolvecollision!(p, o, T::Function, θ::Function, new_ω::Function)
This is the ray-splitting implementation. The three functions given are drawn from
the ray-splitting dictionary that is passed directly to `evolve!()`. For a calculated
incidence angle φ, if T(φ) > rand(), ray-splitting occurs.
"""
function resolvecollision!(p::AbstractParticle, o::Obstacle)::Void
  dist = distance(p, o)

  if dist < 0.0
    relocate!(p, o, dist)
  end
  # Perform specular reflection:
  specular!(p, o)
  return
end

function resolvecollision!(p::AbstractParticle, o::PeriodicWall)
  dist = distance(p, o)

  if dist > 0.0
    relocate!(p, o, dist)
  end
  #perform periodicity
  periodicity!(p, o)
  return
end

####################################################
## Straight Propagation
####################################################

"""
```julia
collisiontime(p::AbstractParticle, o::Obstacle)
```
Calculate the collision time between given particle and obstacle.

The funtion chooses the appropriate method depending on the type of particle
(magnetic or not) as well as the type of the obstacle. Returns the time that the
particle, given its current position and Type, must be propagated to reach the
collision point. This time can be given directly to `propagate!(p, time)` which
brings the particle to the collision point.

In the case of magnetic propagation, there are always two possible collisions.
The function internally decides which of the two will occur first, based on the
sign of the angular velocity of the magnetic particle.
"""
function collisiontime(p::Particle, w::Wall)
  n = normalvec(w, p.pos)
  denom = dot(p.vel, n)
  if denom >= 0; return Inf; end
  t = dot(w.sp-p.pos, n)/denom
  t <= 0.0 ? Inf : t
end

function collisiontime(p::Particle, d::Circular)

  dotp = dot(p.vel, normalvec(d, p.pos))
  # Gotta rethink thins for ray spliting
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

function collisiontime(p::Particle, d::Antidot)

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
    #this makes the phi error come back:
    # if t <= 1e-12
    #   t = -B + sqrtD
    # end
  else
    t = -B + sqrtD
  end

  # If collision time is negative, return Inf:
  t <= 0.0 ? Inf : t
end

"""
```julia
propagate!(p::AbstractParticle, t)
```
Propagate the particle `p` for given time `t`, changing appropriately the the
`p.pos` and `p.vel` fields.

For a `Particle` the propagation is a straight line
(i.e. velocity vector is constant). For a `MagneticParticle` the propagation
is circular motion with cyclic frequency `p.omega` and radius `1/p.omega`.
"""
function propagate!(p::Particle, t::Real)
  # Set initial conditions
  vx0=p.vel[1]
  vy0=p.vel[2]
  # Set current (final) values for `pos` (`vel` does not change)
  p.pos += [vx0*t, vy0*t]
end

"""
    evolve!(p::AbstractParticle, bt::Vector{Obstacle}, t [, ray_splitter])
Evolve the given particle `p` inside the billiard table `bt`. If `t` is of type
`Float64`, evolve for as much time as `t`. If however `t` is of type `Int`,
evolve for as many collisions as `t`.

Return the states of the particle between collisions.

The evolution takes into account the particle's Type.
E.g. if `typeof(p) == MagneticParticle` then magnetic evolution will take place.

## Calling
Call the function like:
```julia
t, poss, vels [, ω]* = evolve!(p, bt, t)
```
(see "Returns" section for more)

To get the position, velocity and time timeseries from the above output,
use the function `construct`:
`xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ttotal)...)`

## Returns
As noted by the "!" at the end of the function, it changes its argument
`p` (particle). Most importantly however, this function also returns the main output
expected by a billiard system. This output is a tuple of 3 (or 4) vectors:
* `ct::Vector{Float64}` : Collision times.
* `poss::Vector{SVector{2}}` : Positions during collisions.
* `vels::Vector{SVector{2}})` : Velocities **exactly after** the collisions.
* `ω`, either `Float64` or `Vector{Float64}` : Angular velocity(/ies).

In the case of straight propagation, only the first 3 are returned.
In the case of magnetic propagation, the 4th value is returned as well.
This is either the angular velocity of the particle (`Float64`), or in the case of
ray-splitting it is a vector of the angular velocities at each time step (`Vector`).

The time `ct[i]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting
from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since
`poss[1], vels[1]` are the initial conditions. The angular velocity `ω[i]` is the one
the particle has while propagating from state `poss[i], vels[i]` to `i+1`.

Notice that at any point, the velocity vector `vels[i]` is the one obtained **after**
the specular reflection of the (i-1)th collision.
The function `construct` takes that into account.

## Ray-splitting billiards
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
function evolve!(p::Particle, bt::Vector{Obstacle}, t)

  if t <= 0
    error("`evolve!()` cannot evolve backwards in time.")
  end

  rt = Float64[]
  rpos = SVector{2,Float64}[]
  rvel = SVector{2,Float64}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, 0.0)

  count = zero(t)

  colobst_idx = 1
  t_to_write = 0.0

  while count < t
    # Declare these because `bt` is of un-stable type!
    tcol::Float64 = 0.0
    tmin::Float64 = Inf

    for i in eachindex(bt)
      tcol = collisiontime(p, bt[i])
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst_idx = i
      end
    end#obstacle loop

    propagate!(p, tmin)
    resolvecollision!(p, bt[colobst_idx])

    t_to_write += tmin

    if typeof(bt[colobst_idx]) <: PeriodicWall
      continue
    else
      push!(rpos, p.pos + p.current_cell)
      push!(rvel, p.vel)
      push!(rt, t_to_write)
      # set counter
      count += increment_counter(t, t_to_write)
      t_to_write = 0.0
    end
  end#time loop
  return (rt, rpos, rvel)
end

"""
```julia
construct(ct, poss, vels[, ω][, dt=0.01])
```
Given the main output of this package (see `evolve!()` function) construct the
timeseries of the position and velocity, as well as the time vector.

In case of not given ω (or ω == 0), straight construction takes place.
In case of ω != 0 or ω::Vector{Real} magnetic construction takes place.

The additional optional argument of `dt` (only valid for Magnetic construction)
is the timestep with which the timeseries are constructed.

## Calling
Call this function like:
```julia
# Straight propagation:
xt, yt, vxt, vyt, ts = construct(ct, poss, vels)
# Magnetic propagation:
xt, yt, vxt, vyt, ts = construct(ct, poss, vels, ω, dt)
# Any-kind-of propagation (skip dt for straight):
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ttotal[, ray_splitter])..., dt)
```
There is no difference in the `construct` call for Ray-Splitting billiards if one
uses the awesome ellipsis operator of Julia. The reason for that is that `evolve!()`
also returns the vector of angular velocities when necessary.

## Returns
A tuple of the following:
* xt::Vector{Float64} : x position time-series
* yt::Vector{Float64} : y position time-series
* vxt::Vector{Float64} : x velocity time-series
* vyt::Vector{Float64} : y velocity time-series
* ts::Vector{Float64} : Continuous time vector
"""
function construct(t::Vector{Float64},
  poss::Vector{SVector{2,Float64}}, vels::Vector{SVector{2,Float64}})

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

function propagate!(p::MagneticParticle, t::Float64)

  # "Initial" conditions
  ω  = p.omega
  vx0= p.vel[1]
  vy0= p.vel[2]
  φ0 = atan2(vy0, vx0)
  # Propagate:
  p.pos += SVector{2, Float64}(sin(ω*t + φ0)/ω - sin(φ0)/ω,
  -cos(ω*t + φ0)/ω + cos(φ0)/ω)
  p.vel = SVector{2, Float64}(cos(ω*t + φ0), sin(ω*t + φ0))
end
propagate!(p::MagneticParticle, t::Real) = propagate!(p, Float64(t))

"""
```julia
realangle(p::MagneticParticle, o::Obstacle, inter::Vector{SVector}, pc, pr)
```
Given the intersections `inter` of the trajectory of a magnetic particle `p` with
some obstacle `o`, find which of the two is the "real" one, i.e. occurs first.
Returns the angle of first collision, which is equal to the time to first
collision divided by ω.

The function also takes care of problems that may arise when particles are very
close to the obstacle's boundaries, due to Floating-point precision.

(the cyclotron center `pc` and radius `pr` are suplimented for efficiency, since they
have been calculated already)
"""
function realangle(p::MagneticParticle, o::Obstacle,
  intersections::Vector{SVector{2, Float64}}, pc::SVector{2, Float64}, pr::Float64)

  ω = p.omega
  P0 = p.pos
  PC = pc - P0
  θ = Inf
  for i in intersections
    d2 = dot(i-P0,i-P0) #distance of particle from intersection point
    # Check dot product for close points:
    if d2 <= 1e-6 #THIS CONDITION MUST DEPEND ON THE Cyclotron Radius!!!
      # AND ALSO ON THE TYPE OF OBSTACLE
      # AND ALSO IF I HAVE TWO INTERSECTIONS ISTEAD OF ONLY ONE!!!
      dotp = dot(p.vel, normalvec(o,  p.pos))
      # Case where velocity points away from obstacle:
      dotp >= 0 && continue
    end

    d2r = (d2/(2pr^2))
    d2r > 2 && (d2r = 2.0)
    # Correct angle value for small arguments (between 0 and π/2):
    θprime = d2r < 1e-8 ? acos1mx(d2r) : acos(1.0 - d2r)

    # Get "side" of i:
    PI = i - P0
    side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
    # Get angle until i (positive number between 0 and 2π)
    side < 0 && (θprime = 2π-θprime)
    # Set minimum angle (first collision)
    if θprime < θ
      θ = θprime
    end
  end
  return θ
end

function collisiontime(p::MagneticParticle, w::Wall)
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
  Δ <= 0.0 && return Inf
  # Intersection coefficients:
  u1 = (-b - sqrt(Δ))/2a
  u2 = (-b + sqrt(Δ))/2a
  cond1 = (0.0 <= u1 <= 1.0)
  cond2 = (0.0 <= u2 <= 1.0)
  # Check if the line is completely inside the circle:
  !cond1 && !cond2 && return Inf
  # Calculate intersection points:
  intersections = SVector{2, Float64}[]
  cond1 && push!(intersections, w.sp + u1*(w.ep - w.sp))
  cond2 && push!(intersections, w.sp + u2*(w.ep - w.sp))
  # Calculate real time until intersection:
  θ = realangle(p, w, intersections, pc, pr)
  # Collision time, equiv. to arc-length until collision point:
  return θ*pr
end

function collisiontime(p::MagneticParticle, o::Circular)
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
  test = rc^2 - a^2
  if test < 0
    println("ω = $ω")
    println("pc = $pc")
    println("rc = $rc")
    println("p1 = $p1")
    println("r1 = $r1")
    println("centers d = $d")
    println("distance = $(distance(p, o))")
    println("a = $a")
    println("rc^2 - a^2 = $test")
  end

  h = sqrt(rc^2 - a^2)
  # Collision points (always 2):
  I1 = SVector{2, Float64}(
  pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
  pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
  I2 = SVector{2, Float64}(
  pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
  pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
  # Calculate real time until intersection:
  θ = realangle(p, o, [I1, I2], pc, rc)
  # Collision time, equiv. to arc-length until collision point:
  return θ*rc
end


function evolve!(p::MagneticParticle, bt::Vector{Obstacle}, t, warning::Bool = false)

  if t <= 0
    error("`evolve!()` cannot evolve backwards in time.")
  end

  ω = p.omega
  absω = abs(ω)
  rt = Float64[]
  rpos = SVector{2,Float64}[]
  rvel = SVector{2,Float64}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, zero(Float64))

  count = zero(t)
  t_to_write = 0.0
  colobst_idx = 1

  while count < t
    tmin::Float64 = Inf
    tcol::Float64 = 0.0

    for i in eachindex(bt)
      tcol = collisiontime(p, bt[i])
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst_idx = i
      end
    end#obstacle loop

    if tmin == Inf
      warning && warn("Pinned particle in evolve! (Inf col t)")
      push!(rpos, rpos[end])
      push!(rvel, rvel[end])
      push!(rt, Inf)
      return (rt, rpos, rvel, ω)
    end

    propagate!(p, tmin)
    resolvecollision!(p, bt[colobst_idx])

    t_to_write += tmin
    # Write output only if the collision was not made with PeriodicWall
    if typeof(bt[colobst_idx]) == PeriodicWall
      # Pinned particle:
      if t_to_write >= 2π/absω
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
      t_to_write = 0.0
    end
  end#time loop
  return (rt, rpos, rvel, ω)
end

function construct(t::Vector{Float64},  poss::Vector{SVector{2,Float64}},
  vels::Vector{SVector{2,Float64}}, ω::Real, dt=0.01)


  if ω == 0
    return construct(t, poss, vels)
  end
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
