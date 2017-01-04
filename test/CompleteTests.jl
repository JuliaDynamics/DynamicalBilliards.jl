if !isdefined(:collisiontime)
  include("Collisions.jl")
  include("Propagation.jl")
  include("PlotBilliards.jl")
end

print("Tests started at: ")
print(Dates.format(now(), "HH:MM:s"), "\n")


check_circle = false
check_straight_sinai = false
check_magnetic_sinai = true
check_straight_sinai_periodic = false
check_magnetic_sinai_periodic = false
check_magnetic_pinned = true

if check_circle
  println("Currently testing if...")# Test:
  println("-Circle billiard works, and randominside() works")
  println("-collisiontime for Circle works")
  println("-resolvecollision for Particle with")
  println(" Circle works and backpropagates")
  println("-evolve!() works and never there is infinite time")
  println("-the particle is always inside the Circle")
  for r in [0.4, 1.5]
    println("...for r = ", r)
    bt = Obstacle{Float64}[Circle([0.5,0.5], r, "Circle")]
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
end

if check_straight_sinai
  println("Currently testing if...")# Test:
  println("-billiard_sinai works, and randominside() works")
  println("-collisiontime for FiniteWall and Disk works")
  println("-resolvecollision for Particle with")
  println(" FiniteWall and Disk works and backpropagates")
  println("-evolve!() works and never there is infinite time")
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
end

if check_magnetic_sinai
  println("Currently testing if...")# Test:
  println("-ωpropagate!() and ωevolve!() work")
  println(" for particle in Sinai billiard")
  println("-back-propagation for magnetic billiards works")
  println("-particle is always inside the billiard")
  println("-ωconstruct works and timeseries is always inside")
  println("-the collision time is always finite")
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
        p = randominside(bt)

        t, poss, vels = ωevolve!(ω, p, bt, tt)
        if t[end] == Inf
          println("Inf. collision time as a result of leakage...?")
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
        xt, yt, vxt, vyt, ts = ωconstruct(ω, t, poss, vels)
        error_level = 1e-15 #this cannot be less than 1e15
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
end


if check_straight_sinai_periodic
  println("Currently testing if...")# Test:
  println("-randominside() works for periodic billiards")
  println("-collisiontime for PeriodicWall works")
  println("-resolvecollision for Particle with")
  println(" PeriodicWall works and forward-propagates")
  println("-evolve!() works and `pos` is always out of Disk")
  println("-minimum collision time is always >= 1-2r")

  for (r, x, y) in [(0.3, 1.5, 1.0), (0.5, 1.4, 2.2), (0.2, 0.8, 2.2)]
    println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai_periodic(r, x, y)
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
end

if check_magnetic_sinai_periodic
  println("Currently testing if...")# Test:
  println("-ωcollisiontime for PeriodicWall works")
  println("-resolvecollision for Particle with")
  println(" PeriodicWall works and forward-propagates")
  println("-ωevolve!() works and `pos` is always out of Disk")
  println("-minimum collision time is always >= 1-2r")
  # Be sure to choose ω where pinned cannot exist
  for (r, x, y) in [(0.3, 1.5, 1.0), (0.5, 1.4, 2.2)]
    for ω in [0.1, 0.5]
      println("...for (ω,r,x,y) = ", (ω, r, x, y))
      bt = billiard_sinai_periodic(r, x, y)
      xmin, ymin, xmax, ymax = cellsize(bt)
      d = bt[5]
      c = d.c
      tt=4000.0
      partnum = 1000
      invalid = 0
      minddist = min(x, y)

      for i in 1:partnum
        p = randominside(bt)
        ts, poss, vels = ωevolve!(ω, p, bt, tt)
        if ts[end] == Inf
          continue
        end


        error_level = 1e-12 #this huge error comes from the modulo operation
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
end

if check_magnetic_pinned
  println("Currently testing if...")# Test:
  println("-there are not any pinned particles for very small ω")
  # Be sure to choose ω where pinned cannot exist
  for (r, x, y) in [(0.4, 1.0, 1.0)]
    for ω in [0.02, 0.04]
      println("...for ω = ", ω)
      bt = billiard_sinai_periodic(r, x, y)
      tt=10000.0
      partnum = 1000
      pinnednum = 0

      for i in 1:partnum
        p = randominside(bt)
        ts, poss, vels = ωevolve!(ω, p, bt, tt)
        if ts[end] == Inf
          pinnednum += 1
          #error("Pinned particle!")
        end
      end#particle loop
      if pinnednum > 1
        error("too many pinned particles")
      end

    end#omega loop
  end#x,y loop
end






















print("Tests ended at ")
println(Dates.format(now(), "HH:MM:s"))
println("without any errors!!")
