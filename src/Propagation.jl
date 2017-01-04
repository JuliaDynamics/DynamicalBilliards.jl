# Collisions.jl must be loaded BEFORE this

####################################################
## Linearized collision times
####################################################
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
function collisiontime{T<:AbstractFloat}(p::AbstractParticle{T}, w::Wall{T})
  n = normalize(w.normal)
  denom = dot(p.vel, n)
  denom == 0 && return convert(T, Inf)
  t = dot(w.sp-p.pos, n)/denom
  t <= zero(T) ? convert(T, Inf) : t
end

function collisiontime{T<:AbstractFloat}(p::AbstractParticle{T}, d::Circular{T})

  infT = convert(T, Inf)
  dc = p.pos - d.c
  B = dot(p.vel, dc)       #pointing towards circle: B < 0
  C = dot(dc, dc) - d.r^2  #being outside of circle: C > 0
  Δ = B^2 - C

  Δ <= zero(T) && return infT
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

  # # Case where sqrt(D)==B but C!=0
  # if sqrtD == B && C != 0.0
  #   #case of C=0 is covered by the above
  #   if C<0
  #     t = (B>0 ? 1e-15 : 2B)
  #   else
  #     t = (B>0 ? infT : 1e-15)
  #   end
  # end

  # BOU HOU HOU WHY DO I HAVE TO DO THIS...
  t <= 1e-14*one(T) ? infT : t
end

function propagate!{T<:AbstractFloat}(p::Particle{T}, t::T)
  # Set initial conditions
  vx0=p.vel[1]
  vy0=p.vel[2] #use vel instead of vxoft because vel changes at collisions
  # Set current (final) values for `pos` (`vel` is unchanged)
  p.pos += [vx0*t, vy0*t]
  #p.pos += SVector{2, T}(vx0*t, vx0*t)
end


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

    tmin -= 1e-15
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

#REWORK THIIIIS
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

function ωpropagate!{T<:AbstractFloat}(ω::T, p::Particle{T}, t::T)

  # "Initial" conditions
  # x0 = p.pos[1]
  # y0 = p.pos[2]
  vx0= p.vel[1]
  vy0= p.vel[2] # use vel instead of vxoft because vel changes at collisions
  φ0 = atan2(vy0, vx0)

  # Set current (final) values for `pos` and `vel`
  # p.pos = SVector{2, T}(
  # sin(ω*t + φ0)/ω + x0 - sin(φ0)/ω, -cos(ω*t + φ0)/ω + y0 + cos(φ0)/ω )
  p.pos += SVector{2, T}( sin(ω*t + φ0)/ω - sin(φ0)/ω, -cos(ω*t + φ0)/ω + cos(φ0)/ω )
  #p.pos = SVector{2, T}([p.xoft[end], p.yoft[end]])
  p.vel = SVector{2, T}(cos(ω*t + φ0), sin(ω*t + φ0))
end


function ωcollisiontime_OLS{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, w::Wall{T})
  pc, pr = cyclotron(ω, p)
  P0 = p.pos
  P2P1 = w.ep - w.sp
  P1P3 = w.sp - pc
  # Solve quadratic:
  a::T = dot(P2P1, P2P1)
  b::T = 2*dot(P2P1, P1P3)
  c::T = dot(P1P3, P1P3) - pr^2
  Δ = b^2 -4*a*c
  # Check if line is completely outside the circle:
  # (at equal it is tangent)
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

  ### Calculate real time until intersection:
  Pa = pc - P0
  Pa3 = SVector{3, T}(Pa..., 0)
  P03 = SVector{3, T}(P0..., 0)
  θ = convert(T, Inf)
  for i in intersections
    d2 = dot(i-P0,i-P0)
    # If i and P0 are identical, skip this point completely
    #d2 == 0 && continue
    #d2 <=1e-22 && continue #at 1e-22 all tests pass
    d2r = (d2/(2pr^2))
    # The number at the following comparison is crucial for resolving corners
    d2r < 1e-22 && continue
    #θprime = (d2r < 1e-14) ? 1e-14 : acos(1 - d2r)
    θprime = acos(1 - d2r)
    # Get "side" of i:
    I3 = SVector{3, T}(i..., 0)
    side = cross((I3-P03), Pa3)[3]*ω
    # Get "side" of i:
    #side = -( (pc[1] - P0[1])*(i[2] - P0[2]) -   (pc[2] - P0[1])*(i[1] - P0[1]) )*ω
    # Get angle until i (positive number between 0 and π)
    side < 0 && (θprime = abs(2π-θprime))
    # Set minimum angle (first collision)
    # notice that 0.1*1e-15 is 1e-16 which is the minimum precision
    # TRY REMOVING SECOND CONDITION
    if θprime < θ && θprime != 0
      θ = θprime
    end
  end
  # Collision time, equiv. to arc-length until collision point:
  return θ*pr
end

#MY IMPLEMENTATION OF CROSS:
function ωcollisiontime{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, w::Wall{T})
  pc, pr = cyclotron(ω, p)
  P0 = p.pos
  P2P1 = w.ep - w.sp
  P1P3 = w.sp - pc
  # Solve quadratic:
  a::T = dot(P2P1, P2P1)
  b::T = 2*dot(P2P1, P1P3)
  c::T = dot(P1P3, P1P3) - pr^2
  Δ = b^2 -4*a*c
  # Check if line is completely outside the circle:
  # (at equal it is tangent)
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

  ### Calculate real time until intersection:
  PC = pc - P0
  θ = convert(T, Inf)
  for i in intersections
    d2 = dot(i-P0,i-P0)
    # If i and P0 are identical, skip this point completely
    #d2 == 0 && continue
    # IS THIS LINE REALLY NECESSARY ???
    d2 <= 1e-22 && continue #at 1e-22 all tests pass
    d2r = (d2/(2pr^2))
    # The number at the following comparison is crucial for resolving corners
    d2r < 1e-16 && continue
    #θprime = (d2r < 1e-14) ? 1e-14 : acos(1 - d2r)
    θprime = acos(1 - d2r)
    # Get "side" of i:
    PI = i - P0
    side = (PI[1]*PC[2] - PI[2]*PC[1])*ω
    # Get "side" of i:
    #side = -( (pc[1] - P0[1])*(i[2] - P0[2]) -   (pc[2] - P0[1])*(i[1] - P0[1]) )*ω
    # Get angle until i (positive number between 0 and π)
    side < 0 && (θprime = abs(2π-θprime))
    # Set minimum angle (first collision)
    # notice that 0.1*1e-15 is 1e-16 which is the minimum precision
    # TRY REMOVING SECOND CONDITION
    if θprime < θ && θprime != 0
      θ = θprime
    end
  end
  # Collision time, equiv. to arc-length until collision point:
  return θ*pr
end



function ωcollisiontime{T<:AbstractFloat}(ω::T, p::AbstractParticle{T}, d::Circular{T})
  pc, rc = cyclotron(ω, p)
  p1 = d.c
  r1 = d.r
  d = norm(p1-pc)
  if (d >= rc + r1) || (d <= abs(rc-r1))
    return convert(T, Inf)
  end
  # Solve quadratic:
  a = (rc^2 - r1^2 + d^2)/2d
  h = sqrt(rc^2 - a^2)
  # Collision points:
  I1 = SVector{2, T}(
  pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
  pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d)
  I2 = SVector{2, T}(
  pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
  pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d)
  ### Calculate real time until intersection:
  Pa = pc - p.pos # Vector from particle to cyclotron center
  P0 = p.pos
  Pa3 = SVector{3, T}(Pa..., 0)
  P03 = SVector{3, T}(P0..., 0)
  θ = convert(T, Inf)
  for I in (I1, I2)
    d2 = dot(I-P0,I-P0)
    #If I and P0 are identical, skip this point completely (Inf time)
    d2 <= 1e-14 && continue
    # Get angle until I (positive number between 0 and π)
    θprime = acos(1 - (d2/(2rc^2)))
    # Get real angle until I:
    I3 = SVector{3, T}(I..., 0)
    side = cross((I3-P03), Pa3)[3]*ω
    side < 0 && (θprime = 2π-θprime)
    # Set minimum angle (first collision), excluding 0 angles
    # notice that 0.1*1e-15 is 1e-16 which is the minimum precision
    # TRY REMOVING SECOND CONDITION
    if θprime < θ && θprime > 1e-15
      θ = θprime
    end
  end
  # Collision time, equiv. to arc-length until collision point:
  return θ*rc
end


function ωevolve!{T<:AbstractFloat}(ω::T, p::AbstractParticle{T},
  bt::Vector{Obstacle{T}}, ttotal::T, dt::T=0.05*one(T))

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
      tcol = ωcollisiontime(ω, p, obst)
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
    tmin -= 1e-15
    ωpropagate!(ω, p, tmin)
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

function ωconstruct{T<:AbstractFloat}(ω::T, t::Vector{T},
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

    # if length(timevec) == 0
    #   em = "Error! length of timevec = 0 in ωconstruct !!!\n"
    #   em*= "t0=$t0, colt=$colt"
    #   error(em)
    # end

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
