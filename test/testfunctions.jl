
function check_raysplitting_omega(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing: check_raysplitting_omega")# Test:
    println("--ray-splitting with Antidot works")
    println("--for straight propagation")
    println("--for magnetic propagation")
    println("--when transmission NOT depending on ω")
    println("--Also checking `isphysical()`")
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
  x = 1.0; y = 1.0
  bt = billiard_rectangle(x, y)
  a = Antidot([0.5, 0.5], 0.3, true)
  push!(bt, a)
  if !isphysical(rayspl)
    error("Given ray-splitter is not physical!")
  end

  for ω in [0.0, 1.0, -0.5, 0.02]
    printinfo && println("...for ω = ", ω)
    error_level = 1e-9

    for i in 1:partnum
      p = randominside(bt, ω)
      xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0, rayspl)...)
      dmaxx = maximum(xt) - x
      dmaxx > error_level && error("xmax - x = $dmaxx")
      dminx = minimum(xt)
      dminx < -error_level && error("xmin = $dminx")
      dmaxy = maximum(yt) - y
      dmaxy > error_level && error("ymax - y = $dmaxy")
      dminy = minimum(yt)
      dminy < -error_level && error("ymin = $dminy")
      bt[5].pflag = true
    end
  end
  return true
end

function check_raysplitting_periodic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing: check_raysplitting_periodic")# Test:
    println("--ray-splitting with Antidot works in periodic billiard")
    println("--for straight propagation")
    println("--for magnetic propagation")
    println("--when transmission does NOT depend on ω")
    println("--Also testing `isphysical()`")
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
      bt[5].pflag = true
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
    println("--for both straight and magnetic propagation")
  end
  if partnum > 10
    partnum = 10
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
    println("--resolvecollision for AbstractParticle with")
    println("  random obstacle works and reflects randomly")
    println("--evolve!() works and never there is infinite time")
  end
  for (r, x, y) in [(0.25, 1.0, 1.0), (0.35, 1.4, 2.2)]
    printinfo && println("...for (r,x,y) = ", (r, x, y))
    bt = billiard_sinai(r, x, y; setting = "random")
    d = bt[5]
    c = d.c
    tt=1000.0
    for ω in [0.0, 0.02, 0.78596]
      printinfo && println("......and ω = ", ω)
      for i in 1:partnum
        p = randominside(bt, ω)
        ts, poss, vels = evolve!(p, bt, tt)
        if ts[end] == Inf
          println("Final pos:")
          println("px = $(p.pos[1])")
          println("py = $(p.pos[2])")
          println("Final vel:")
          println("vx = $(p.vel[1])")
          println("vy = $(p.vel[2])")
          println("Length of ts = $(length(ts))")
          # plot_billiard(bt)
          # plot_particle(p)
          # return p, bt
          error("Infinite collision time in Random Sinai")
        end

        error_level = 1e-15

        xt = [pos[1] for pos in poss]; yt = [pos[2] for pos in poss]
        dist = sqrt.(((xt .- c[1]).^2 .+ (yt .- c[2]).^2))
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
    end#omega loop
  end#x,y loop
  return true
end



function check_klein_magnetic(partnum; printinfo = true)
  if printinfo
    println("\nCurrently testing: check_klein_magnetic")# Test:
    println("--ray-splitting with Antidot works in a periodic billiard")
    println("--For complicated, magnetic field dependent Tunneling")
    println("--Emulates klein tunneling in magnetic fields")
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
      ct, ps, vs = evolve!(p, bt, 1000.0, rayspl)
      if ct[end] == Inf
        error("Infinite collision time in periodic sinai with Antidot (pinned)!")
      end
    end
  end
  return true
end


function check_lyapunov_spectrum(partnum; printinfo = true)
    if printinfo
        println("\nCurrently testing: check_lyapunov_spectrum")
    println("--lyapunovspectrum works for a polygonar lorentz  gas with periodic walls")
  end

    l = 2.0
    r = 1.0
    sides = 6
    printinfo && println("...for n = $(sides) sides, l_polygon = $l and r_disk = $r  ")
    bt =  bt = billiard_polygon(6,l; setting = "periodic")
    disc = Disk([0., 0.], r)
    push!(bt, disc)
    tt=10000.0

    for i in 1:partnum
        p = randominside(bt)
        exps = lyapunovspectrum(p, bt, tt)
        error_level = 1e-2
        sumpair = exps[1] + exps[4]

        abs(sumpair) > error_level && error("λ1 + λ4 = $(sumpair) > $error_level")

        exps[2] > error_level && error("λ2 > $error_level")

        exps[3] > error_level && error("λ3 > $error_level")
    end#particle loop
    return true
end
