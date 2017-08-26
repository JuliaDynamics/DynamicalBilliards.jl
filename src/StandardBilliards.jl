export billiard_rectangle, billiard_sinai, billiard_polygon, billiard_lorentz

####################################################
## Famous/Standard Billiards
####################################################
"""
```julia
billiard_rectangle(x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard table allow for ray-splitting.
"""
function billiard_rectangle(x=1.0, y=1.0; setting::String = "standard")

  bt = Obstacle[]
  if setting == "standard"
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
  elseif setting == "periodic"
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
  elseif setting == "random"
    o = 0.0
    sp = [o,o]; ep = [o, y]; n = [x,o]
    leftw = RandomWall(sp, ep, n, "Left random wall")
    sp = [x,o]; ep = [x, y]; n = [-x,o]
    rightw = RandomWall(sp, ep, n, "Right random wall")
    sp = [o,y]; ep = [x, y]; n = [o,-y]
    topw = RandomWall(sp, ep, n, "Top random wall")
    sp = [o,o]; ep = [x, o]; n = [o,y]
    botw = RandomWall(sp, ep, n, "Bottom random wall")
    push!(bt, leftw, rightw, topw, botw)
  elseif setting == "ray-splitting"
    o = 0.0
    sp = [o,o]; ep = [o, y]; n = [x,o]
    leftw = SplitterWall(sp, ep, n, "Left ray-splitting wall")
    sp = [x,o]; ep = [x, y]; n = [-x,o]
    rightw = SplitterWall(sp, ep, n, "Right ray-splitting wall")
    sp = [o,y]; ep = [x, y]; n = [o,-y]
    topw = SplitterWall(sp, ep, n, "Top ray-splitting wall")
    sp = [o,o]; ep = [x, o]; n = [o,y]
    botw = SplitterWall(sp, ep, n, "Bottom ray-splitting wall")
    push!(bt, leftw, rightw, topw, botw)
  else
    error("The given setting=$setting is unknown.")
  end
  return bt
end

"""
```julia
billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
```
Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type,
  enforcing periodicity at the boundaries
* "random" : The velocity is randomized upon collision.
* "ray-splitting" : All obstacles in the billiard table allow for ray-splitting.
"""
function billiard_sinai(r=0.25, x=1.0, y=1.0; setting = "standard")
  if (setting == "periodic") && (r>=x/2 || r>=y/2)
    es = "Disk radius too big for a periodic Sinai billiard.\n"
    es*= "Obstacles must not overlap with `PeriodicWall`s."
    error(es)
  end
  bt = billiard_rectangle(x, y; setting = setting)
  c = [x/2, y/2]
  if setting == "random"
    centerdisk = RandomDisk(c, r, "Random disk")
  elseif setting == "ray-splitting"
    centerdisk = Antidot(c, r, "Antidot")
  else
    centerdisk = Disk(c, r, "Disk")
  end
  push!(bt, centerdisk)
end

"""
    billiard_lorentz(r=0.25, x=1.0, y=1.0)
Alias for `billiard_sinai(r,x,y; setting = "periodic")`.
"""
billiard_lorentz(r=0.25, x=1.0, y=1.0) = billiard_sinai(r,x,y; setting = "periodic")

"""
    billiard_polygon(n::Int, R, center = [0,0]; setting = "standard")
Return a vector of obstacles that defines a regular-polygonal billiard table
with `n` sides, radius `r` and given `center`.

Note: `R` denotes the so-called outer radius, not the inner one.

### Settings
* "standard" : Specular reflection occurs during collision.
* "periodic" : The walls are `PeriodicWall` type, enforcing periodicity
  at the boundaries. Only available for `n=4` or `n=6`.
* "random" : The velocity is randomized upon collision.
"""
function billiard_polygon(sides::Int, r::Real, center = [0,0]; setting = "standard")
  bt = Obstacle[]
  verteces = [[r*cos(2π*i/sides), r*sin(2π*i/sides)] + center for i in 1:sides]

  if setting == "standard"
    T = FiniteWall
    wallname = "wall"
  elseif setting == "periodic"
    if sides != 4 && sides != 6
      error("Polygonal and periodic billiard can exist only for `n=4` or `n=6`")
    end
    T = PeriodicWall
    wallname = "periodic wall"
    inr = r*cos(π/sides)
  elseif setting == "random"
    T = RandomWall
    wallname = "random wall"
  end

  for i in eachindex(verteces)
    N = length(verteces)
    starting = verteces[i]
    ending = verteces[mod1(i+1, N)]
    # Normal vector must look at where the particle is coming from
    w = ending - starting
    if setting == "periodic"
      normal = inr*normalize([-w[2], w[1]])
    else
      normal = [-w[2], w[1]]
    end
    wall = FiniteWall(starting, ending, normal, wallname*" $i")
    push!(bt, wall)
  end
  return bt
end
