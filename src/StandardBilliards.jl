export billiard_rectangle, billiard_sinai, billiard_polygon

####################################################
## Famous/Standard Billiards
####################################################
"""
    billiard_rectangle(x=1.0, y=1.0; periodic = false)
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).
"""
function billiard_rectangle(x=1.0, y=1.0; periodic = false)

  bt = Obstacle[]
  if !periodic
    o = 0.0
    sp = [o,o]; ep = [o, y]; n = [x,o]
    leftw2 = FiniteWall(sp, ep, n, "Left wall")
    sp = [x,o]; ep = [x, y]; n = [-x,o]
    rightw2 = FiniteWall(sp, ep, n, "Right wall")
    sp = [o,y]; ep = [x, y]; n = [o,-y]
    topw2 = FiniteWall(sp, ep, n, "Top wall")
    sp = [o,o]; ep = [x, o]; n = [o,y]
    botw2 = FiniteWall(sp, ep, n, "Bottom wall")
    push!(bt, leftw2, rightw2, topw2, botw2)
  else
    o = 0.0
    sp = [o,o]; ep = [o, y]; n = [x,o]
    leftw = PeriodicWall(sp, ep, n, "Left periodic boundary")
    sp = [x,o]; ep = [x, y]; n = [-x,o]
    rightw = PeriodicWall(sp, ep, n, "Right periodic boundary")
    sp = [o,y]; ep = [x, y]; n = [o,-y]
    topw = PeriodicWall(sp, ep, n, "Top periodic boundary")
    sp = [o,o]; ep = [x, o]; n = [o,y]
    botw = PeriodicWall(sp, ep, n, "Bottom periodic boundary")
    push!(bt, leftw, rightw, topw, botw)
  end
  return bt
end

"""
    billiard_sinai(r, x=1.0, y=1.0; periodic = false)
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".
"""
function billiard_sinai(r, x=1.0, y=1.0; periodic = false)
  if (periodic == true) && (r>=x/2 || r>=y/2)
    es = "Disk radius too big for a periodic Sinai billiard.\n"
    es*= "Obstacles must not overlap with `PeriodicWall`s."
    error(es)
  end
  bt = billiard_rectangle(x,y; periodic = periodic)
  c = [x/2, y/2]
  centerdisk = Disk(c, r, "Disk")
  push!(bt, centerdisk)
end

"""
    billiard_polygon(n::Int, R, center = [0,0]; periodic = true)
Return a vector of obstacles that defines a regular-polygonal billiard table
with `n` sides, radius `r` and given `center`. If `n` is even, you may choose a
periodic version of the billiard.

Note: `R` denotes the so-called outer radius, not the inner one.
"""
function billiard_polygon(sides::Int, r::Real, center = [0,0]; periodic = false)
  bt = Obstacle[]
  verteces = [[r*cos(2π*i/sides), r*sin(2π*i/sides)] + center for i in 1:sides]
  if !periodic
    for i in eachindex(verteces)
      N = length(verteces)
      starting = verteces[i]
      ending = verteces[mod1(i+1, N)]
      # Normal vector must look at where the particle is coming from
      w = ending - starting
      normal = [-w[2], w[1]]
      wall = FiniteWall(starting, ending, normal, "Wall $i")
      push!(bt, wall)
    end
  else
    !iseven(sides) && error("A periodic billiard must have even number of sides.")
    for i in eachindex(verteces)
      N = length(verteces)
      starting = verteces[i]
      ending = verteces[mod1(i+1, N)]
      # Normal vector must look at where the particle is coming from
      # and in the case of periodic must have length exactly as much as it is
      # from one side to the opposite
      w = ending - starting
      inr = r*cos(π/sides)
      normal = inr*normalize([-w[2], w[1]])
      wall = PeriodicWall(starting, ending, normal, "Periodic wall $i")
      push!(bt, wall)
    end
  end
  return bt
end
