if !isdefined(:collisiontime)
  include("ParticlesObstacles.jl")
  include("Propagation.jl")
  include("StandardBilliards.jl")
  include("PlotBilliards.jl")
  include("RaySplitting.jl")
end

print("Tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")

function check_circle()
  println("Currently testing if...")# Test:
  println("--Circle billiard works, and randominside() works")
  println("--collisiontime for Circle works")
  println("--resolvecollision for Particle with")
  println("  Circle works and backpropagates")
  println("--evolve!() works and never there is infinite time")
  println("--the particle is always inside the Circle")
  for r in [0.4, 1.5]
    println("...for r = ", r)
    bt = Obstacle[a = Antidot([0.5, 0.5], 0.5, false, "circle")]
    d = bt[1]
    c = d.c
    tt=10000.0
    partnum = 1000

    for i in 1:partnum
      p = randominside(bt)
      ts, poss, vels = evolve!(p, bt, tt)

      error_level = 1e-16

      xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
      dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
      mind = maximum(dist)
      if mind - d.r > error_level
        error("min((x,y)-r) = $(mind - d.r)")
      end

    end#particle loop
  end#x,y loop
  return true
end

function check_straight_sinai()
  println("Currently testing if...")# Test:
  println("--billiard_sinai works, and randominside() works")
  println("--collisiontime for FiniteWall and Disk works")
  println("--resolvecollision for Particle with")
  println("  FiniteWall and Disk works and backpropagates")
  println("--evolve!() works and never there is infinite time")
  for (r, x, y) in [(0.4, 1.0, 1.0), (0.3, 1.5, 1.0), (0.5, 1.4, 2.2)]
    println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y)
    d = bt[5]
    c = d.c
    tt=10000.0
    partnum = 1000

    for i in 1:partnum
      p = randominside(bt)
      ts, poss, vels = evolve!(p, bt, tt)

      error_level = 1e-16

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

function check_magnetic_sinai()
  println("Currently testing if...")# Test:
  println("--ωpropagate!() and ωevolve!() work")
  println("  for particle in Sinai billiard")
  println("--back-propagation for magnetic billiards works")
  println("--particle is always inside the billiard")
  println("--ωconstruct works and timeseries is always inside")
  println("--the collision time is always finite")
  for ω in [-0.5, 1.2]
    for (r, x, y) in [(0.4, 1.0, 1.0), (0.3, 1.5, 1.0)]
      print("...for ω = $ω and ")
      print("for (r,x,y) = ", (r, x, y), "\n")
      bt = billiard_sinai(r, x, y)
      d = bt[5]
      c = d.c
      tt=1000.0
      partnum = 1000
      for i in 1:partnum
        p = randominside(ω, bt)

        t, poss, vels = evolve!(p, bt, tt)
        if t[end] == Inf
          println("Inf. collision time as a result of leakage or corner.")
          println("p.pos = $(p.pos)")
          error("Pinned particle in closed sinai billiard... (inf col. time)")
        end

        xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
        error_level = 0.0
        dist = sqrt(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
        mind = minimum(dist)
        if mind - d.r < error_level
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


function check_straight_sinai_periodic()
  println("Currently testing if...")# Test:
  println("--randominside() works for periodic billiards")
  println("--collisiontime for PeriodicWall works")
  println("--resolvecollision for Particle with")
  println("  PeriodicWall works and forward-propagates")
  println("--evolve!() works and `pos` is always out of Disk")
  println("--minimum collision time is always >= 1-2r")

  for (r, x, y) in [(0.3, 1.5, 1.0), (0.5, 1.4, 2.2), (0.2, 0.8, 2.2)]
    println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y; periodic=true)
    xmin, ymin, xmax, ymax = cellsize(bt)
    d = bt[5]
    c = d.c
    tt=10000.0
    partnum = 1000
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

function check_magnetic_sinai_periodic()
  println("Currently testing if...")# Test:
  println("--ωcollisiontime for PeriodicWall works")
  println("--resolvecollision for Particle with")
  println("  PeriodicWall works and forward-propagates")
  println("--ωevolve!() works and `pos` is always out of Disk")
  println("--minimum collision time is always >= 1-2r")
  # Be sure to choose ω where pinned cannot exist
  for (r, x, y) in [(0.4, 1.5, 1.0), (0.5, 1.4, 2.2)]
    for ω in [0.1, 0.5]
      println("...for (ω,r,x,y) = ", (ω, r, x, y))
      bt = billiard_sinai(r, x, y; periodic=true)
      xmin, ymin, xmax, ymax = cellsize(bt)
      d = bt[5]
      c = d.c
      tt=4000.0
      partnum = 1000
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

function check_magnetic_pinned()
  println("Currently testing if...")# Test:
  println("--there are not any pinned particles for very small ω")
  # Be sure to choose ω where pinned cannot exist
  for (r, x, y) in [(0.4, 1.0, 1.0)]
    for ω in [0.02, 0.04]
      println("...for (ω,r,x,y) = ", (ω, r, x, y))
      bt = billiard_sinai(r, x, y; periodic = true)
      tt=10000.0
      partnum = 2000

      for i in 1:partnum
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

function check_previous_obstacle()
  println("Currently testing if...")# Test:
  println("--The previous collision obstacle is never the same as the")
  println("  current in straight closed sinai billiard.")
  println("--Minimum collision time is never less than 1e-10")
  ttotal = 10000.0
  bt = billiard_sinai(0.3)

  for j in 1:1000
    p = randominside(bt)
    colobst = bt[1]
    prev_obst = nothing
    tcount = 0.0
    t_to_write = 0.0

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
        println("Current obstacle: $(colobst.name)")
            println("tmin = $tmin")
        error("Previus obstacle same as current obstacle")
      end
      if tcount !=0 && tmin < 1e-10
        error("Minimum collision time = $tmin")
      end

      propagate!(p, tmin)
      dt = resolvecollision!(p, colobst)
      t_to_write += tmin + dt
      tcount += t_to_write
      t_to_write = 0.0
      prev_obst = colobst
    end#time loop
  end#particle loop
  return true
end#test



function check_raysplitting_omega()
  println("Currently testing if...")# Test:
  println("- ray-splitting with antidot works")
  println("  for magnetic and straight propagation")

  sa = (θ, where, ω) -> where ? 2.0*θ : 0.5*θ
  T = (θ, where, ω) -> begin
    if where
      abs(θ) < π/4 ? 0.5*exp(-(θ)^2/2(π/8)^2) : 0.0
    else
      0.5*exp(-(θ)^2/2(π/4)^2)
    end
  end
  newo = ((x, bool) -> bool ? -0.5x : 2.0x)
  rayspl = Dict{Int, Vector{Function}}(5 => [T, sa, newo])

  bt = billiard_rectangle()
  a = Antidot([0.5, 0.5], 0.3, true)
  push!(bt, a)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 1.0, -0.5, 0.02]
    println("...for ω = ", ω)


    for i in 1:500
      p = randominside(bt, ω)
      xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0, rayspl)...)
    end
  end
  return true
end

function check_raysplitting_periodic()
  println("Currently testing if...")# Test:
  println("- ray-splitting with antidot works in periodic billiard")
  println("  for magnetic and straight propagation")

  sa = (θ, where, ω) -> where ? 2.0*θ : 0.5*θ
  T = (θ, where, ω) -> begin
    if where
      abs(θ) < π/4 ? 0.5*exp(-(θ)^2/2(π/8)^2) : 0.0
    else
      0.5*exp(-(θ)^2/2(π/4)^2)
    end
  end
  newo = ((x, bool) -> bool ? -0.5x : 2.0x)
  rayspl = Dict{Int, Vector{Function}}(5 => [T, sa, newo])

  bt = billiard_rectangle(periodic = true)
  a = Antidot([0.5, 0.5], 0.3, true)
  push!(bt, a)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 1.0, -0.5, 0.02]
    println("...for ω = ", ω)


    for i in 1:500
      p = randominside(bt, ω)
      xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0, rayspl)...)
    end
  end
  return true
end



#check_circle()
check_straight_sinai()
check_magnetic_sinai()
check_straight_sinai_periodic()
check_magnetic_sinai_periodic()
check_magnetic_pinned()
check_previous_obstacle()
check_raysplitting_omega()
check_raysplitting_periodic()















print("Tests ended at ")
println(Dates.format(now(), "HH:MM:s"))
println("without any errors!!")
