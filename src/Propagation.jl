# Collisions.jl must be loaded BEFORE this

function acos1mx(x)
  sqrt(2x) + sqrt(x)^3/(6sqrt(2))
end

####################################################
## Resolve Collisions
####################################################
#All these functions will become dependend on velocity angle for ray-splitting billiards.
"""
    specular!(p::AbstractParticle, o::Obstacle)
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
function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T},
  printstuff = false)

  #println("Resolving collision with $(o.name)")
  #id = randstring(4)

  dist = distance(p, o)
  dt = zero(T)
  #printstuff = false
  if dist < 0.0
    # if abs(dist) > 1e-8
    #   println("In resolve with $(o.name), dist = $dist")
    #   error("Distance unreasonably big")
    # end
    if printstuff
      println("In resolve with $(o.name), dist = $dist")
    end
    #dt = 10lct(p, o, dist) #linearized collision time, should be negative
    dt = 10lct(p, o, dist) #linearized collision time, should be negative
    printstuff && println("2*l.c.t. = $dt")
    if dt > 0
      n = normalvec(o, p.pos)
      println("dot(p.vel, n) = $(dot(p.vel, n))")
      println("(should be negative)")
      error("lct should be negative")
    end
    # Check the orders of magnitude:
    vxt=p.vel[1]*dt; vyt= p.vel[2]*dt
    #exponent of velocity*timeL
    ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
    #exponent of position
    pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
    # Ensure that the position will be changed
    if pe - ve > 16.0
      dt *= 10^(pe - ve - 16.0)
    end
    #println("dt = $dt\n")
    # Propagate backwards:
    p.pos += [p.vel[1]*dt, p.vel[2]*dt]
  elseif printstuff
    println("Good distance for $(o.name), dist = $dist")
  end

  #TEST TO BE SURE
  if distance(p, o)<0.0
    error("Dist = $(distance(p,o)) after back-propagation")
  end


  # Perform specular reflection:
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
  printstuff && println("-----------------------------------------------------\n")
  return dt
end

function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::PeriodicWall{T},
  printstuff = false)

    #println("Resolving collision with $(o.name)")
    #id = randstring(4)
    dist = distance(p, o)

    t = zero(T)
    #printstuff = false

    if dist > 0.0
      # if abs(dist) > 1e-8
      #   println("In resolve with $(o.name), dist = $dist")
      #   error("Distance unreasonably big")
      # end
      if printstuff

        println("In resolve with $(o.name), dist = $dist")
      end
      t = 10lct(p, o, dist) #linearized collision time, should be positive
      if t < 0
        n = normalvec(o, p.pos)
        println("dot(p.vel, n) = $(dot(p.vel, n))")
        println("(should be negative)")
        error("lct should be positive")
      end
      # Check the orders of magnitude:
      vxt=p.vel[1]*t; vyt= p.vel[2]*t
      #exponent of velocity*timeL
      ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
      #exponent of position
      pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
      # Ensure that the position will be changed
      if pe - ve > 16
        t *= 10^(pe - ve - 16)
      end
      # Propagate forwards:
      p.pos += [p.vel[1]*t, p.vel[2]*t]
    elseif printstuff
      println("Good distance for $(o.name), dist = $dist")
    end
  #TEST TO BE SURE
  if distance(p, o)>0.0
    error("Dist = $(distance(p,o)) after forward-propagation")
  end

  #perform periodicity
  p.pos += o.normal
  p.current_cell -= o.normal
  printstuff && println("-----------------------------------------------------\n")
  return t
end


####################################################
## Linearized collision times
####################################################
"""
    lct(p::AbstractParticle, o::Obstacle, distance)
Return the Linearized Collision Time (` (+/-) dist/dot(p.vel, normal)`) between particle
and obstacle, given the calculated distance between them.
"""
function lct{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T}, dist::T)
  n = normalvec(o, p.pos)
  t = -dist/dot(p.vel, n)
end


function lct{T<:AbstractFloat}(p::AbstractParticle{T}, w::Circle{T}, dist::T)
  n = normalvec(w, p.pos)
  t = dist/dot(p.vel, n)
end

####################################################
## Straight Propagation
####################################################

"""
    collisiontime(p::AbstractParticle, o::Obstacle)
Calculate the (minimum) collision time between given particle and obstacle.

The funtion chooses the appropriate method depending on the type of particle (magnetic or
not) as well as the type of the obstacle. Returns the time that the particle, given its
current position, must be propagated to reach the collision point.

In the case of magnetic propagation, there are always two possible collisions. The function
internally decides which of the two will occur first, based on the sign of the angular
velocity of the magnetic particle.
"""
function collisiontime{T<:AbstractFloat}(p::Particle{T}, w::Wall{T})
  n = normalize(w.normal)
  denom = dot(p.vel, n)
  denom == 0 && return convert(T, Inf)
  t = dot(w.sp-p.pos, n)/denom
  t <= zero(T) ? convert(T, Inf) : t
end

function collisiontime{T<:AbstractFloat}(p::Particle{T}, d::Circular{T})

  infT = convert(T, Inf)
  dc = p.pos - d.c
  B = dot(p.vel, dc)       #pointing towards circle: B < 0
  C = dot(dc, dc) - d.r^2  #being outside of circle: C > 0
  Δ = B^2 - C

  Δ <= 0 && return infT
  sqrtD = sqrt(Δ)

  # Case of being slightly outside and looking outside:
  # if B > 0.0 && C > 0.0
  #   return infT
  # end
  # Closest point:
  t = -B - sqrtD
  # Case of being on top of looking inside:
  if t==0.0 && B < 0.0
    t=-2B
  # Case of being inside but closest point is in negative time
  elseif t < 0.0 && C < 0.0
    t = -B + sqrtD
  end

  # BOU HOU HOU WHY DO I HAVE TO DO THIS...
  t <= 1e-14*one(T) ? infT : t
end

"""
    propagate!(p::AbstractParticle, t)
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
    (t::Vector, poss::Vector{SVector{2}}, vels:: Vector{SVector{2}})
As noted by the "!" at the end of the function, it changes its argument `p` (particle).
Most importantly however, this function also returns the main output expected by a billiard
system. This output is a typle of three vectors:
* `t::Vector` : Collision times.
* `poss::Vector{SVector{2}}` : Positions during collisions.
* `vels:: Vector{SVector{2}})` : Velocities **exactly after** the collisions.

The time `t[i]` is the time necessary to reach state `poss[i], vels[i]` starting from the
state `poss[i-1], vels[i-1]`. That is why `t[0]` is always 0 since `poss[0], vels[0]` are
the initial conditions.

Notice that at any point, the velocity vector `vels[i]` is the one obtained **after**
the specular reflection of the (i-1)th collision.
The function `construct` takes that into account.

# Usage in Ray-splitting billiards

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
      if tcol < 1e-10
        println("Collision time with $(obst.name):")
        println(" -> tcol = $tcol")
      end
      #DELETTHIS
      if tcol == 0.0 #This is important!!!
        es = "tcol was returned 0.0 in evolve!\n"
        es*= "calculated for col, with $(obst.name)\n"
        es*= "prev. col. was made with $(prev_obst.name)"
        error(es)
      end
      # println("Collision time with $(obst.name):")
      # println(" -> tcol = $tcol")
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst = obst
      end
    end#obstacle loop

    if tmin == infT
      error("Collision time infinite; Impossible error in evolve!")
    end

    # tmin -= 1e-15
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
  #here you call check_validity(p, bt) if you want to,
  #but resolve_collision does it anyway
  # isin = distance(p, bt)
  #
  # if isin == -1.0 || isin == 0.0
  #   em = "AFTER evolve, distance with $(colobst.name): $isin \n"
  #   em *="p.pos = $(p.pos)\n"
  #   em *="x[end] = $(p.xoft[end]), "
  #   em *="y[end] = $(p.yoft[end])\n"
  #   em *="It last collided with $(colobst.name)\n"
  #   error(em)
  # end
  return (rt, rpos, rvel)
end

"""
```
construct{T<:AbstractFloat}(t::Vector{T},
          poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}}, dt=0.5*one(T))

construct(omega, t, poss, vels, dt=0.01*one(T))
```
Given the main output of this package (see `evolve!()` function) construct the timeseries
of the position and velocity, as well as the time vector.

# Calling
Call this function like:
`xt, yt, vxt, vyt, ts = construct(t, poss, vels)`
for straight propagation, and:
`xt, yt, vxt, vyt, ts = construct(omega, t, poss, vels)`
for magnetic propagation, with `omega = p.omega` the angular velocity of the particle.

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
  vxt= [vel[1] for vel in vels]
  vyt= [vel[2] for vel in vels]
  # CREATE CORRECT T! ! !
  return xt, yt, vxt, vyt, t
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
```
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
    # MAYBE 1e-10 is too small???? and in the corners it will lead to problems?
    if d2 <= 1e-10
      dotp = dot(p.vel, normalvec(o,  p.pos))
      # Case where velocity points away from obstacle:
      dotp > 0 && continue
    end

    d2r = (d2/(2pr^2))
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
  ### Calculate real time until intersection:
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
  push!(rt, zero(T)) #using 0.0 instead works perfectly fine

  infT = convert(T, Inf)
  colobst = NullWall(T)
  tcount = zero(T)
  t_to_write = zero(T)

  while tcount <= ttotal
    tmin = infT

    for obst in bt
      tcol = collisiontime(p, obst)
      # println("Collision time with $(obst.name):")
      # println(" -> tcol = $tcol")
      #DELETE THIS !!!
      if tcol == 0.0 #This is important!!! But only for ωevolve!?
        es = "tcol was returned 0.0 in ωevolve!\n"
        es*= "calculated for col, with $(obst.name)\n"
        error(es)
      end
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

    # Makin gthe colision time a bit shorter reduces many computations
    # Because makes almost all distances "good"
    # tmin -= 1e-15
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

    if colt < 0
      error("colt < 0, major error, destruction.
      Probably it was backpropagated too much when the collision time was also too small")
    end

    t0 = ct[i-1]
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
## Magnetic Propagation OLD STUFF
####################################################
#
# function ωpropagate!{T<:AbstractFloat}(ω::T, p::Particle{T}, t::T)
#
#   # "Initial" conditions
#   # x0 = p.pos[1]
#   # y0 = p.pos[2]
#   vx0= p.vel[1]
#   vy0= p.vel[2] # use vel instead of vxoft because vel changes at collisions
#   φ0 = atan2(vy0, vx0)
#
#   # Set current (final) values for `pos` and `vel`
#   # p.pos = SVector{2, T}(
#   # sin(ω*t + φ0)/ω + x0 - sin(φ0)/ω, -cos(ω*t + φ0)/ω + y0 + cos(φ0)/ω )
#   p.pos += SVector{2, T}( sin(ω*t + φ0)/ω - sin(φ0)/ω, -cos(ω*t + φ0)/ω + cos(φ0)/ω )
#   #p.pos = SVector{2, T}([p.xoft[end], p.yoft[end]])
#   p.vel = SVector{2, T}(cos(ω*t + φ0), sin(ω*t + φ0))
# end
#
#
# function realangle{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T},
#   intersections::Vector{SVector{2, T}}, ω::T, pc::SVector{2, T}, pr::T)
#
#   P0 = p.pos
#   PC = pc - P0
#   θ = convert(T, Inf)
#   for i in intersections
#     d2 = dot(i-P0,i-P0)
#
#     # Check dot product for close points:
#     # MAYBE 1e-10 is too small???? and in the corners it will lead to problems?
#     if d2 <= 1e-10
#       dotp = dot(p.vel, normalvec(o,  p.pos))
#       # Case where velocity points away from obstacle:
#       dotp > 0 && continue
#     end
#
#     d2r = (d2/(2pr^2))
#     θprime = d2r < 1e-16 ? acos1mx(d2r) : acos(1-d2r)
#
#     # Get "side" of i:
#     PI = i - P0
#     side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
#     # Get angle until i (positive number between 0 and π)
#     side < 0 && (θprime = abs(2π-θprime))
#     # Set minimum angle (first collision)
#     if θprime < θ
#       θ = θprime
#     end
#   end
#   return θ
# end
#
#
# # IMPLEMENT CHECKING Dot of velocity with normal:
# function ωcollisiontime{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, w::Wall{T})
#   pc, pr = cyclotron(ω, p)
#   P0 = p.pos
#   P2P1 = w.ep - w.sp
#   P1P3 = w.sp - pc
#   # Solve quadratic:
#   a::T = dot(P2P1, P2P1)
#   b::T = 2*dot(P2P1, P1P3)
#   c::T = dot(P1P3, P1P3) - pr^2
#   Δ = b^2 -4*a*c
#   # Check if line is completely outside the circle:
#   # (at equal it is tangent)
#   Δ <= zero(T) && return convert(T, Inf)
#   # Intersection coefficients:
#   u1 = (-b - sqrt(Δ))/2a
#   u2 = (-b + sqrt(Δ))/2a
#   cond1 = (zero(T) <= u1 <= one(T))
#   cond2 = (zero(T) <= u2 <= one(T))
#   # Check if the line is completely inside the circle:
#   !cond1 && !cond2 && return convert(T, Inf)
#   # Calculate intersection points:
#   intersections = SVector{2, T}[]
#   cond1 && push!(intersections, w.sp + u1*(w.ep - w.sp))
#   cond2 && push!(intersections, w.sp + u2*(w.ep - w.sp))
#
#   ### Calculate real time until intersection:
# #  PC = pc - P0
# #  θ = convert(T, Inf)
#   θ = realangle(p, w, intersections, ω, pc, pr)
#   # for i in intersections
#   #   d2 = dot(i-P0,i-P0)
#   #
#   #   # Check dot product for close points:
#   #   # MAYBE 1e-10 is too small???? and in the corners it will lead to problems?
#   #   if d2 <= 1e-10
#   #     dotp = dot(p.vel, normalvec(w,  p.pos))
#   #     # Case where velocity points away from obstacle:
#   #     dotp > 0 && continue
#   #   end
#   #
#   #   d2r = (d2/(2pr^2))
#   #   if d2r < 1e-16
#   #     θprime = acos1mx(d2r)
#   #   else
#   #     θprime = acos(1-d2r)
#   #   end
#   #
#   #   # Get "side" of i:
#   #   PI = i - P0
#   #   side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
#   #   # Get angle until i (positive number between 0 and π)
#   #   side < 0 && (θprime = abs(2π-θprime))
#   #   # Set minimum angle (first collision)
#   #   if θprime < θ
#   #     θ = θprime
#   #   end
#   # end
#   # Collision time, equiv. to arc-length until collision point:
#   return θ*pr
# end
#
#
# function ωcollisiontime_OLD{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, d::Circular{T})
#   pc, rc = cyclotron(ω, p)
#   p1 = d.c
#   r1 = d.r
#   d = norm(p1-pc)
#   if (d >= rc + r1) || (d <= abs(rc-r1))
#     return convert(T, Inf)
#   end
#   # Solve quadratic:
#   a = (rc^2 - r1^2 + d^2)/2d
#   h = sqrt(rc^2 - a^2)
#   # Collision points:
#   I1 = SVector{2, T}(
#   pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
#   pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
#   I2 = SVector{2, T}(
#   pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
#   pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
#   ### Calculate real time until intersection:
#   PC = pc - p.pos
#   P0 = p.pos
#   θ = convert(T, Inf)
#   for i in (I1, I2)
#     d2 = dot(i-P0,i-P0)
#     #If I and P0 are identical, skip this point completely (Inf time)
#     #d2 <= 1e-14 && continue
#     d2r = (d2/(2rc^2))
#     # The number at the following comparison is crucial for resolving corners
#     d2r < 1e-16 && continue
#     # Get angle until I (positive number between 0 and π)
#     θprime = acos(1 - d2r)
#     # Get real angle until I:
#     PI = i - P0
#     side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
#     side < 0 && (θprime = 2π-θprime)
#     # Set minimum angle (first collision), excluding 0 angles
#     # notice that 0.1*1e-15 is 1e-16 which is the minimum precision
#     # TRY REMOVING SECOND CONDITION
#     if θprime < θ && θprime > 1e-8
#       θ = θprime
#     end
#   end
#   # Collision time, equiv. to arc-length until collision point:
#   return θ*rc
# end
#
# function ωcollisiontime{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, o::Circular{T})
#   pc, rc = cyclotron(ω, p)
#   p1 = o.c
#   r1 = o.r
#   d = norm(p1-pc)
#   if (d >= rc + r1) || (d <= abs(rc-r1))
#     return convert(T, Inf)
#   end
#   # Solve quadratic:
#   a = (rc^2 - r1^2 + d^2)/2d
#   h = sqrt(rc^2 - a^2)
#   # Collision points:
#   I1 = SVector{2, T}(
#   pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
#   pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
#   I2 = SVector{2, T}(
#   pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
#   pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
#   ### Calculate real time until intersection:
#   θ = realangle(p, o, [I1, I2], ω, pc, rc)
#   # Collision time, equiv. to arc-length until collision point:
#   return θ*rc
# end
#
#
# function ωevolve!{T<:AbstractFloat}(ω::T, p::AbstractParticle{T},
#   bt::Vector{Obstacle{T}}, ttotal::T, dt::T=0.05*one(T))
#
#   rt = T[]
#   rpos = SVector{2,T}[]
#   rvel = SVector{2,T}[]
#   push!(rpos, p.pos)
#   push!(rvel, p.vel)
#   push!(rt, zero(T)) #using 0.0 instead works perfectly fine
#
#   infT = convert(T, Inf)
#   colobst = NullWall(T)
#   tcount = zero(T)
#   t_to_write = zero(T)
#
#   while tcount <= ttotal
#     tmin = infT
#
#     for obst in bt
#       tcol = ωcollisiontime(ω, p, obst)
#       # println("Collision time with $(obst.name):")
#       # println(" -> tcol = $tcol")
#       #DELETE THIS !!!
#       if tcol == 0.0 #This is important!!! But only for ωevolve!?
#         es = "tcol was returned 0.0 in ωevolve!\n"
#         es*= "calculated for col, with $(obst.name)\n"
#         error(es)
#       end
#       # Set minimum time:
#       if tcol < tmin
#         tmin = tcol
#         colobst = obst
#       end
#     end#obstacle loop
#
#     if tmin == infT
#       println("pinned particle! (Inf col t)")
#       push!(rpos, rpos[end])
#       push!(rvel, rvel[end])
#       push!(rt, convert(T, Inf))
#       return (rt, rpos, rvel)
#     end
#
#     # Makin gthe colision time a bit shorter reduces many computations
#     # Because makes almost all distances "good"
#     # tmin -= 1e-15
#     ωpropagate!(ω, p, tmin)
#     dt = resolvecollision!(p, colobst)
#     t_to_write += tmin + dt
#     # Write output only if the collision was not made with PeriodicWall
#     if typeof(colobst) <: PeriodicWall
#       # Pinned particle:
#       if t_to_write >= 2π/ω
#         println("pinned particle! (completed circle)")
#         push!(rpos, rpos[end])
#         push!(rvel, rvel[end])
#         push!(rt, convert(T, Inf))
#         return (rt, rpos, rvel)
#       end
#       #If not pinned, continue (do not write for PeriodicWall)
#       continue
#     else
#       push!(rpos, p.pos + p.current_cell)
#       push!(rvel, p.vel)
#       push!(rt, t_to_write)
#       tcount += t_to_write
#       t_to_write = zero(T)
#     end
#
#   end#time loop
#   return (rt, rpos, rvel)
# end
#
# function ωconstruct{T<:AbstractFloat}(ω::T, t::Vector{T},
#   poss::Vector{SVector{2,T}}, vels::Vector{SVector{2,T}}, dt=0.01*one(T))
#
#   xt = [poss[1][1]]
#   yt = [poss[1][2]]
#   vxt= [vels[1][1]]
#   vyt= [vels[1][2]]
#   ts = [t[1]]
#   ct = cumsum(t)
#
#   for i in 2:length(t)
#     φ0 = atan2(vels[i-1][2], vels[i-1][1])
#     x0 = poss[i-1][1]; y0 = poss[i-1][2]
#     colt=t[i]
#
#     if colt < 0
#       error("colt < 0, major error, destruction.
#       Probably it was backpropagated too much when the collision time was also too small")
#     end
#
#     t0 = ct[i-1]
#     if colt >= dt
#       timevec = collect(0:dt:colt)[2:end]
#       timevec[end] == colt || push!(timevec, colt)
#     else
#       timevec = colt
#     end
#
#     # if length(timevec) == 0
#     #   em = "Error! length of timevec = 0 in ωconstruct !!!\n"
#     #   em*= "t0=$t0, colt=$colt"
#     #   error(em)
#     # end
#
#     for td in timevec
#       push!(vxt, cos(ω*td + φ0))
#       push!(vyt, sin(ω*td + φ0))
#       push!(xt, sin(ω*td + φ0)/ω + x0 - sin(φ0)/ω)  #vy0 is sin(φ0)
#       push!(yt, -cos(ω*td + φ0)/ω + y0 + cos(φ0)/ω) #vx0 is cos(φ0)
#       push!(ts, t0 + td)
#     end#collision time
#   end#total time
#   return xt, yt, vxt, vyt, ts
# end
