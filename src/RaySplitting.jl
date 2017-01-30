export isphysical, acceptable_raysplitter, supports_raysplitting

"""
```julia
supports_raysplitting(obst::Obstacle)
```
Return `true` if the given obstacle supports ray-splitting.
"""
function supports_raysplitting(obst::Obstacle)
  n = fieldnames(typeof(obst))
  in(:where, n) ? true : false
end

"""
```julia
acceptable_raysplitter(raysplitter, bt)
```
Check if the given ray-splitting dictionary `raysplitter` can be used in conjuction with
given billiard table `bt`.
"""
function acceptable_raysplitter(ray::Dict{Int, Vector{Function}}, bt::Vector{Obstacle})
  for i in keys(ray)
    if !supports_raysplitting(bt[i])
      println("Obstacle at index $i of given billiard table does not have a field `where`")
      println("and therefore does not support ray-splitting.")
      return false
    end
  end
  true
end


# Resolve collision for Ray-splitting
function resolvecollision!(p::MagneticParticle, a::Obstacle, T::Function,
  θ::Function, new_ω::Function = ((x, bool) -> x))

  dt = 0.0
  ω = p.omega
  # Determine incidence angle (0 < θ < π/4)
  n = normalvec(a, p.pos)
  φ = acos(dot(p.vel, -n))
  # if this is wrong then my normal vec is wrong:
  # if φ >= π/2
  #   println("φ=$φ")
  #   println("Collision happens with $(a.name)")
  #   println("Distance is: $(distance(p, a))")
  #   #error("φ shoud be between 0 and π/2")
  #   println("-----------")
  # end
  # ray-splitting (step 2)
  if T(φ, a.where, ω) > rand()
    # Step 3
    if cross2D(p.vel, n) < 0
      φ *= -1
    end
    # Step 4
    theta = θ(φ, a.where, ω)
    # Step 5
    a.where = !a.where
    # Step 6
    n = normalvec(a, p.pos) #notice that this is reversed! It's the new!
    Θ = theta + atan2(n[2], n[1])
    # Step 7
    dist = distance(p, a)  #this is also reversed! It's the new distance!
    # Step 8
    if dist < 0.0
      dt = relocate!(p, a, dist)
    end
    # Step 9
    p.vel = [cos(Θ), sin(Θ)]
    # Step 10
    p.omega = new_ω(ω, !a.where)  # notice the exclamation mark !
  # No ray-splitting:
  else
    dist = distance(p, a)
    dt = 0.0

    if dist < 0.0
      dt = relocate!(p, a, dist)
    end
    #perform specular
    specular!(p, a)
  end
  if abs(dt) > 1e-4
    error("dt = $dt in resolve ray-splitting Magnetic.")
  end
  return dt
end

# Resolve collision for ray-splitting
function resolvecollision!(p::Particle, a::Obstacle, T::Function, θ::Function)

  dt = 0.0
  ω = 0.0
  # Determine incidence angle (0 < θ < π/4)
  n = normalvec(a, p.pos)
  φ = acos(dot(p.vel, -n))
  # if this is wrong then my normal vec is wrong:
  if φ >= π/2
    println("φ=$φ")
    if a.where == true
      println("Particle should be coming from outside of disk")
    else
      println("Particle should be coming from inside of disk")
    end
    error("φ shoud be between 0 and π/2")
  end
  # ray-splitting (step 2)
  if T(φ, a.where, ω) > rand()
    # Step 3
    if cross2D(p.vel, n) < 0
      φ *= -1
    end
    # Step 4
    theta = θ(φ, a.where, ω)
    # Step 5
    a.where = !a.where
    # Step 6
    n = normalvec(a, p.pos) #notice that this is reversed! It's the new!
    Θ = theta + atan2(n[2], n[1])
    # Step 7
    dist = distance(p, a)  #this is also reversed! It's the new distance!
    # Step 8
    if dist < 0.0
      dt = relocate!(p, a, dist)
    end
    # Step 9
    p.vel = [cos(Θ), sin(Θ)]
  # No ray-splitting:
  else
    dist = distance(p, a)
    if dist < 0.0
      dt = relocate!(p, a, dist)
    end
    #perform specular
    specular!(p, a)
  end
  if abs(dt) > 1e-4
    error("dt = $dt in resolve ray-splitting Straight.")
  end
  return dt
end

# For Particle and Ray-Splitting:
function evolve!(p::Particle, bt::Vector{Obstacle}, ttotal::Float64,
  ray::Dict{Int, Vector{Function}})

  rt = Float64[]
  rpos = SVector{2,Float64}[]
  rvel = SVector{2,Float64}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, 0.0)
  tcount = 0.0
  colobst = bt[1]
  prev_obst = bt[1]
  colind::Int = length(bt)
  t_to_write = 0.0


  while tcount < ttotal
    tmin = Inf

    for i in eachindex(bt)
      obst = bt[i]
      tcol = collisiontime(p, obst)
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst = obst
        colind = i
      end
    end#obstacle loop
    # if tmin < 1e-10 && tcount!=0
    #   println("-----------------")
    #   println("In raysplit evolve, tmin = $tmin")
    #   println("tcount = $tcount")
    #   println("Collision is to happen with $(colobst.name)")
    #   println("Antidot state: where = $(bt[5].where)")
    #   println("Previous collision with: $(prev_obst.name)")
    #
    # end

    propagate!(p, tmin)
    if haskey(ray, colind)
      dt = resolvecollision!(p, colobst, ray[colind][1], ray[colind][2])
    else
      dt = resolvecollision!(p, colobst)
    end
    t_to_write += tmin + dt
    prev_obst = colobst

    if typeof(colobst) == PeriodicWall
      continue
    else
      push!(rpos, p.pos + p.current_cell)
      push!(rvel, p.vel)
      push!(rt, t_to_write)
      tcount += t_to_write
      t_to_write = 0.0
    end
  end#time loop
  return (rt, rpos, rvel)
end

# For MagneticParticle and Ray-Splitting. Returns one extra vector with omegas!!!
function evolve!(p::MagneticParticle, bt::Vector{Obstacle},
                 ttotal::Float64, ray::Dict{Int, Vector{Function}})

  omegas = Float64[]
  rt = Float64[]
  rpos = SVector{2,Float64}[]
  rvel = SVector{2,Float64}[]
  push!(rpos, p.pos)
  push!(rvel, p.vel)
  push!(rt, zero(Float64))
  push!(omegas, p.omega)

  tcount = 0.0
  t_to_write = 0.0
  colobst = bt[1]
  prev_obst = bt[end]
  colind = 1

  while tcount < ttotal
    tmin = Inf

    for i in eachindex(bt)
      obst = bt[i]
      tcol = collisiontime(p, obst)
      # Set minimum time:
      if tcol < tmin
        tmin = tcol
        colobst = obst
        colind = i
      end
    end#obstacle loop

    if tmin == Inf
      println("pinned particle! (Inf col t)")
      push!(rpos, rpos[end])
      push!(rvel, rvel[end])
      push!(rt, Inf)
      return (rt, rpos, rvel)
    end
    if tmin < 1e-10 && tcount!=0
      println("-----------------")
      println("In raysplit evolve, tmin = $tmin")
      println("tcount = $tcount")
      println("Collision is to happen with $(colobst.name)")
      println("$(colobst.name) state: where = $(colobst.where)")
      println("Previous collision with: $(prev_obst.name)")
    end

    propagate!(p, tmin)
    if haskey(ray, colind)
      dt = resolvecollision!(p, colobst, ray[colind][1], ray[colind][2], ray[colind][3])
    else
      dt = resolvecollision!(p, colobst)
    end
    t_to_write += tmin + dt
    # Write output only if the collision was not made with PeriodicWall
    if typeof(colobst) == PeriodicWall
      # Pinned particle:
      if t_to_write >= 2π/abs(p.omega)
        #println("t_to_write = $t_to_write while circle time = $")
        println("pinned particle! (completed circle)")
        push!(rpos, rpos[end])
        push!(rvel, rvel[end])
        push!(rt, Inf)
        push!(omegas, p.omega)
        return (rt, rpos, rvel, omegas)
      end
      #If not pinned, continue (do not write for PeriodicWall)
      prev_obst = colobst
      continue
    else
      push!(rpos, p.pos + p.current_cell)
      push!(rvel, p.vel)
      push!(rt, t_to_write)
      push!(omegas, p.omega)
      tcount += t_to_write
      t_to_write = 0.0
    end
    # if (prev_obst == colobst) && (typeof(colobst) == SplitterWall)
    #   error("Double collision with splitter wall!\n tmin was $tmin\ndt was $dt")
    # end
    #prev_obst = colobst
  end#time loop
  return (rt, rpos, rvel, omegas)
end

function construct(t::Vector{Float64}, poss::Vector{SVector{2,Float64}},
  vels::Vector{SVector{2,Float64}}, omegas::Vector{Float64}, dt=0.01)

  xt = [poss[1][1]]
  yt = [poss[1][2]]
  vxt= [vels[1][1]]
  vyt= [vels[1][2]]
  ts = [t[1]]
  ct = cumsum(t)

  for i in 2:length(t)
    ω = omegas[i-1]
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

"""
    isphysical(raysplitter::Dict{Int, Vector{Function}}; only_mandatory = false)
Return `true` if the given ray-splitting dictionary represends the physical world.
Specifically, check if (φ is the incidence angle):
* Critical angle means total reflection: If θ(φ) ≥ π/2 then T(φ) = 0
* Transmission probability is even function: T(φ) ≈ T(-φ)
* Refraction angle is odd function: θ(φ) ≈ -θ(-φ)
* Ray reversal is true: θ(θ(φ, where, ω), !where, ω) ≈ φ
* Magnetic conservation is true: (ω_new(ω_new(ω, where), !where) ≈ ω
The first property is mandatory and must hold for correct propagation.
They keyword `only_mandatory` notes whether the rest of
the properties should be tested or not.
"""
function isphysical(ray::Dict{Int, Vector{Function}}; only_mandatory = false)
  for i in keys(ray)
    scatter = ray[i][2]
    tr = ray[i][1]
    om = ray[i][3]
    range = -π/2:0.001:π/2
    orange = -1.0:0.01:1.0
    display_er = true
    for where in [true, false]
      for ω in orange
        for φ in range
          θ::Float64 = 0.0
          # Calculate refraction angle:
          try
            θ = scatter(φ, where, ω)
          catch er
            if display_er
              println("Error message: $er")
              println("while calculating the refraction angle with settings:")
              println("index = $i, φ = $φ, where = $where, ω = $ω")
            end
            display_er = false
            T = tr(φ, where, ω)
            if T!= 0
              println("Error message: $er")
              println("while calculating the refraction angle with settings:")
              println("index = $i, φ = $φ, where = $where, ω = $ω")
              println("Also, transmission prob. was not 0 for these settings. (major problem!)")
              return false
            else
              continue
            end
          end
          # Calculate transmission probability:
          T = tr(φ, where, ω)
          # Check critical angle:
          if θ >= π/2 && T > 0
            es = "Refraction angle >= π/2 and T > 0 !\n"
            es*= "For index = $i, tested with φ = $φ, where = $where, ω = $ω"
            println(es)
            return false
          end
          if !only_mandatory
            # Check symmetry:
            if !isapprox(θ, -scatter(-φ, where, ω))
              es = "Scattering angle function is not odd!\n"
              es *="For index = $i, tested with φ = $φ, where = $where, ω = $ω"
              println(es)
              return false
            end
            if !isapprox(T, tr(-φ, where, ω))
              es = "Transmission probability function is not even!\n"
              es *="For index = $i, tested with φ = $φ, where = $where, ω = $ω"
              println(es)
              return false
            end
            # Check ray-reversal:
            if !isapprox(scatter(θ, !where, ω), φ)
              es = "Ray-reversal does not hold!\n"
              es *="For index = $i, tested with φ = $φ, where = $where, ω = $ω"
              println(es)
              return false
            end
            if !isapprox(om(om(ω, where), !where), ω)
              es = "Magnetic reversal does not hold!\n"
              es *="For index = $i, tested with φ = $φ, where = $where, ω = $ω"
              println(es)
              return false
            end
          end
        end#φ range
      end#ω range
    end#where range
  end#obstacle range
  return true
end
