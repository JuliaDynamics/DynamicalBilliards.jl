function check_straight_sinai(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--billiard_sinai works, and randominside() works")
    println("--collisiontime for FiniteWall and Disk works")
    println("--resolvecollision for Particle with")
    println("  FiniteWall and Disk works and backpropagates")
    println("--evolve!() works and never there is infinite time")
  end
  for (r, x, y) in [(0.4, 1.0, 1.0), (0.3, 1.5, 1.0), (0.5, 1.4, 2.2)]
    printinfo && println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y)
    d = bt[5]
    c = d.c
    tt=10000.0

    for i in 1:partnum
      p = randominside(bt)
      ts, poss, vels = evolve!(p, bt, tt)

      error_level = 1e-15

      xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
      dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
      mind = minimum(dist)
      if mind - d.r < -error_level
        error("min((x,y)-r) = $(mind - d.r)")
      end
      dmaxx = maximum(xt) - x
      dmaxx > error_level && error("xmax - x = $dmaxx")
      dminx = minimum(xt)
      dminx < -error_level && error("xmin = $dminx")
      dmaxy = maximum(yt) - y
      dmaxy > error_level && error("ymax - y = $dmaxy")
      dminy = minimum(yt)
      dminy < -error_level && error("ymin = $dminy")
    end#particle loop
  end#x,y loop
  return true
end

function check_magnetic_sinai(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ωpropagate!() and ωevolve!() work")
    println("  for particle in Sinai billiard")
    println("--back-propagation for magnetic billiards works")
    println("--particle is always inside the billiard")
    println("--ωconstruct works and timeseries is always inside")
    println("--the collision time is always finite")
  end
  for ω in [-0.5, 1.2]
    for (r, x, y) in [(0.4, 1.0, 1.0), (0.3, 1.5, 1.0)]
      printinfo && print("...for ω = $ω and ")
      printinfo && print("for (r,x,y) = ", (r, x, y), "\n")
      bt = billiard_sinai(r, x, y)
      d = bt[5]
      c = d.c
      tt=1000.0
      for i in 1:partnum
        p = randominside(ω, bt)

        t, poss, vels = evolve!(p, bt, tt)
        if t[end] == Inf
          println("Inf. collision time as a result of leakage or corner.")
          println("p.pos = $(p.pos)")
          error("Pinned particle in closed sinai billiard... (inf col. time)")
        end

        xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
        error_level = 1e-15
        dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
        mind = minimum(dist)
        if mind - d.r < -error_level
          error("POS: min((x,y)-r) = $(mind - d.r)")
        end
        dmaxx = maximum(xt) - x
        dmaxx > error_level && error("POS: xmax - x = $dmaxx")
        dminx = minimum(xt)
        dminx < -error_level && error("POS: xmin = $dminx")
        dmaxy = maximum(yt) - y
        dmaxy > error_level && error("POS: ymax - y = $dmaxy")
        dminy = minimum(yt)
        dminy < -error_level && error("POS: ymin = $dminy")

        # TWO ERROR CHECKS. One with pos and level = 0.0
        # one with xt and error level 15
        xt, yt, vxt, vyt, ts = construct(t, poss, vels, ω)
        #this cannot be less than 1e15 due to resolvecoolision!() and construct()
        error_level = 1e-15
        dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
        mind = minimum(dist)
        if mind - d.r < -error_level
          error("min((x,y)-r) = $(mind - d.r)")
        end
        dmaxx = maximum(xt) - x
        dmaxx > error_level && error("xmax - x = $dmaxx")
        dminx = minimum(xt)
        dminx < -error_level && error("xmin = $dminx")
        dmaxy = maximum(yt) - y
        dmaxy > error_level && error("ymax - y = $dmaxy")
        dminy = minimum(yt)
        dminy < -error_level && error("ymin = $dminy")
      end#particle loop
    end#radius loop
  end#omega loop
  return true
end


function check_straight_sinai_periodic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--randominside() works for periodic billiards")
    println("--collisiontime for PeriodicWall works")
    println("--resolvecollision for Particle with")
    println("  PeriodicWall works and forward-propagates")
    println("--evolve!() works and `pos` is always out of Disk")
    println("--minimum collision time is always >= 1-2r")
  end

  for (r, x, y) in [(0.3, 1.5, 1.0), (0.5, 1.4, 2.2), (0.2, 0.8, 2.2)]
    printinfo && println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y; setting="periodic")
    xmin, ymin, xmax, ymax = cellsize(bt)
    d = bt[5]
    c = d.c
    tt=10000.0
    invalid = 0
    minddist = min(x, y)

    for i in 1:partnum
      mincolt = Inf
      p = randominside(bt)
      ts, poss, vels = evolve!(p, bt, tt)

      if length(ts) > 2
        mincolt = minimum(ts[3:end])
      else
        continue
      end

      error_level = 1e-8 #this huge error comes from the modulo operation
      xt = [mod(pos[1], xmax) for pos in poss]
      yt = [mod(pos[2], ymax) for pos in poss]

      dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
      mind = minimum(dist)
      if mind - d.r < -error_level
        error("min((x,y)-r) = $(mind - d.r)")
      end

      if mincolt < minddist-2r
        println("min. col. t = $mincolt")
        error("Min. col. time less than disk-disk distance (= $(minddist-2r)")
      end

    end#particle loop
  end#x,y loop
  return true
end

function check_magnetic_sinai_periodic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ωcollisiontime for PeriodicWall works")
    println("--resolvecollision for Particle with")
    println("  PeriodicWall works and forward-propagates")
    println("--ωevolve!() works and `pos` is always out of Disk")
    println("--minimum collision time is always >= 1-2r")
  end
  # Pinned particles are excluded, so don't care about ω
  for (r, x, y) in [(0.4, 1.5, 1.0), (0.25, 1.0, 1.0)]
    for ω in [0.1, 2.18]
      printinfo && println("...for (ω,r,x,y) = ", (ω, r, x, y))
      bt = billiard_sinai(r, x, y; setting="periodic")
      xmin, ymin, xmax, ymax = cellsize(bt)
      d = bt[5]
      c = d.c
      tt=4000.0
      invalid = 0
      minddist = min(x, y)

      for i in 1:partnum
        p = randominside(ω, bt)
        ts, poss, vels = evolve!(p, bt, tt)
        if ts[end] == Inf
          continue
        end


        error_level = 1e-10 #this huge error comes from the modulo operation
        xt = [mod(pos[1], xmax) for pos in poss]
        yt = [mod(pos[2], ymax) for pos in poss]

        dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
        mind = minimum(dist)
        if mind - d.r < -error_level
          error("min((x,y)-r) = $(mind - d.r)")
        end
        mincolt = minimum(ts[3:end])
        if mincolt < minddist-2r
          println("min. col. t = $mincolt")
          error("Min. col. time less than disk-disk distance (= $(minddist-2r)")
        end

      end#particle loop
    end#omega loop
  end#x,y loop
  return true
end

function check_magnetic_pinned(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--there are not any pinned particles for very small ω")
  end
  # Be sure to choose ω pflag pinned cannot exist
  for (r, x, y) in [(0.4, 1.0, 1.0)]
    for ω in [0.02, 0.04]
      printinfo && println("...for (ω,r,x,y) = ", (ω, r, x, y))
      bt = billiard_sinai(r, x, y; setting="periodic")
      tt=1000.0

      for i in 1:2partnum
        p = randominside(ω, bt)
        ts, poss, vels = evolve!(p, bt, tt)
        if ts[end] == Inf
          error("Pinned particle for ω=$ω ")
        end
      end#particle loop
    end#omega loop
  end#x,y loop
  return true
end

function check_previous_obstacle(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--The previous collision obstacle is never the same as the")
    println("  current in a closed sinai billiard.")
    println("--For straight propagation")
    println("--For magnetic propagation with very small fields")
    println("  (same obstacle not possible for very small distances)")
  end
  ttotal = 1000.0
  bt = billiard_sinai(0.3)

  for ω in [0.0, 0.002, -0.004]
    printinfo && println("...for ω=$ω")

    for j in 1:partnum
      p = randominside(bt, ω)
      d = distance(p, bt) #initial distance
      prev_pos = p.pos
      prev_vel = p.vel
      new_pos = p.pos
      new_vel = p.vel
      colobst = nothing
      prev_obst = nothing
      tcount = 0.0
      t_to_write = 0.0
      colnumber = 0

      while tcount < ttotal
        tmin = Inf

        for obst in bt
          tcol = collisiontime(p, obst)
          # Set minimum time:
          if tcol < tmin
            tmin = tcol
            colobst = obst
          end
        end#obstacle loop

        if colobst == prev_obst
          if ω==0
            println("Collision with obstacle: $(colobst.name)")
            println("collision time = $tmin")
            println("init. distance = $d")
            println("current pos = $(new_pos)")
            println("current vel = $(new_vel)")
            println("previus pos = $(prev_pos)")
            println("previus vel = $(prev_vel)")
            println("Collision number: $colnumber")
            error("Previuus obstacle was same as current for straight prop.")
          else
            if tmin<1e-6
              println("Collision with obstacle: $(colobst.name)")
              println("collision time = $tmin")
              println("init. distance = $d")
              println("current pos = $(new_pos)")
              println("current vel = $(new_vel)")
              println("previus pos = $(prev_pos)")
              println("previus vel = $(prev_vel)")
              println("Collision number: $colnumber")
              error("Previuus obstacle was same as current for straight prop.")
            end
          end
        end

        propagate!(p, tmin)
        dt = resolvecollision!(p, colobst)
        t_to_write += tmin + dt
        tcount += t_to_write
        t_to_write = 0.0
        prev_obst = colobst
        prev_pos = new_pos
        prev_vel = new_vel
        new_pos=p.pos
        new_vel=p.vel
        colnumber += 1


      end#time loop
    end#particle loop
  end#omega loop
  return true
end#test



function check_raysplitting_omega(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ray-splitting with Antidot works")
    println("  for magnetic and straight propagation")
  end

  sa = (θ, pflag, ω) -> pflag ? 2.0*θ : 0.5*θ
  T = (θ, pflag, ω) -> begin
    if pflag
      abs(θ) < π/4 ? 0.5*exp(-(θ)^2/2(π/8)^2) : 0.0
    else
      0.5*exp(-(θ)^2/2(π/4)^2)
    end
  end
  newo = ((x, bool) -> bool ? -0.5x : -2.0x)
  rayspl = Dict{Int, Vector{Function}}(5 => [T, sa, newo])

  bt = billiard_rectangle()
  a = Antidot([0.5, 0.5], 0.3, true)
  push!(bt, a)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 1.0, -0.5, 0.02]
    printinfo && println("...for ω = ", ω)


    for i in 1:partnum
      p = randominside(bt, ω)
      xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0, rayspl)...)
    end
  end
  return true
end

function check_raysplitting_periodic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ray-splitting with Antidot works in periodic billiard")
    println("  for magnetic and straight propagation")
  end

  sa = (θ, pflag, ω) -> pflag ? 2.0*θ : 0.5*θ
  T = (θ, pflag, ω) -> begin
    if pflag
      abs(θ) < π/4 ? 0.5*exp(-(θ)^2/2(π/8)^2) : 0.0
    else
      0.5*exp(-(θ)^2/2(π/4)^2)
    end
  end
  newo = ((x, bool) -> bool ? -2.0x : -0.5x)
  rayspl = Dict{Int, Vector{Function}}(5 => [T, sa, newo])

  bt = billiard_rectangle(setting="periodic")
  a = Antidot([0.5, 0.5], 0.3, true)
  push!(bt, a)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 0.02, -0.04]
    printinfo && println("...for ω = ", ω)


    for i in 1:partnum
      p = randominside(bt, ω)
      if ω == 0
        ct, ps, vs = evolve!(p, bt, 1000.0, rayspl)
      else
        ct, ps, vs, os = evolve!(p, bt, 1000.0, rayspl)
      end
      if ct[end] == Inf
        error("Infinite collision time in periodic sinai with Antidot (pinned)!")
      end
    end
  end
  return true
end

function check_splitterwall(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ray-splitting with SplitterWall and Antidot works")
  end

  sa = (θ, pflag, ω) -> pflag ? 2.0*θ : 0.5*θ
  Tp = (p) -> (θ, pflag, ω) -> begin
    if pflag
      abs(θ) < π/4 ? p*exp(-(θ)^2/2(π/8)^2) : 0.0
    else
      (1-p)*exp(-(θ)^2/2(π/4)^2)
    end
  end
  newo = ((x, bool) -> bool ? -0.5x : -2.0x)
  rayspl = Dict{Int, Vector{Function}}(
  5 => [Tp(0.5), sa, newo],
  6 => [Tp(0.35), sa, newo],
  7 => [Tp(0.65), sa, newo])

  bt = billiard_rectangle(2, 1)
  sw = SplitterWall([1.0, 0.0], [1,1], [-1,0], true)
  push!(bt, sw)
  a1 = Antidot([0.5, 0.5], 0.3, "Left Antidot")
  push!(bt, a1)
  a2 = Antidot([1.5, 0.5], 0.2, "Right Antidot")
  push!(bt, a2)

  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, -0.04, 0.5]
    printinfo && println("...for ω = ", ω)

    for i in 1:partnum
      p = randominside(bt, ω)
      xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0, rayspl)...)
    end
  end
  return true
end

function check_random_sinai(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--billiard_sinai works with setting = `random` ")
    println("--collisiontime for RandomWall and RandomDisk works")
    println("--resolvecollision for Particle with")
    println("  random obstacle works and reflects randomly")
    println("--evolve!() works and never there is infinite time")
  end
  for (r, x, y) in [(0.4, 1.0, 1.0), (0.5, 1.4, 2.2)]
    printinfo && println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y; setting = "random")
    d = bt[5]
    c = d.c
    tt=1000.0

    for i in 1:partnum
      p = randominside(bt)
      ts, poss, vels = evolve!(p, bt, tt)

      error_level = 1e-15

      xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
      dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
      mind = minimum(dist)
      if mind - d.r < -error_level
        error("min((x,y)-r) = $(mind - d.r)")
      end
      dmaxx = maximum(xt) - x
      dmaxx > error_level && error("xmax - x = $dmaxx")
      dminx = minimum(xt)
      dminx < -error_level && error("xmin = $dminx")
      dmaxy = maximum(yt) - y
      dmaxy > error_level && error("ymax - y = $dmaxy")
      dminy = minimum(yt)
      dminy < -error_level && error("ymin = $dminy")
    end#particle loop
  end#x,y loop
  return true
end



function check_klein_magnetic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing if...")# Test:
    println("--ray-splitting with Antidot works in a periodic billiard")
    println("  for very complicated, magnetic field dependent Tunneling")
    println("--emulates klein tunneling in magnetic fields")
  end
  if partnum < 1000
    partnum = 1000
  end

  #Create raysplitting
  α = 300.0;  w = 100.0;  n = 1.0
  B0 = 233.3*sqrt(n)/α #this value is in tesla
  Bstar = 116.6*sqrt(n)/w #this value is in tesla
  Bstar /= B0
  kf = sqrt(π*n)
  a = -π*kf*w

  function Transmission(φ, pflag, ω)
    B = ω/2
    β = - B/Bstar
    γ = 1/sqrt(1 - β^2);

    if pflag == true
      return exp( a*γ^3*(sin(φ) - B/Bstar)^2 )
    else
      return exp( a*γ^3*(sin(φ) + B/Bstar)^2 )
    end
  end

  sangle(φ, pflag, ω) = -φ
  newo(x, bool) = -x

  rayspl = Dict{Int, Vector{Function}}(5 => [Transmission, sangle, newo])

  bt = billiard_rectangle()
  ad = Antidot([0.5, 0.5], 0.25, true)
  push!(bt, ad)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 0.16, -0.5, -1.0]
    printinfo && println("...for ω = ", ω)


    for i in 1:partnum
      p = randominside(bt, ω)
      ct, ps, vs = evolve!(p, bt, 4000.0, rayspl)
      if ct[end] == Inf
        error("Infinite collision time in periodic sinai with Antidot (pinned)!")
      end
    end
  end
  return true
end
