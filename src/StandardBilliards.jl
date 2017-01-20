export billiard_rectangle, billiard_sinai

####################################################
## Famous/Standard Billiards
####################################################
"""
    billiard_rectangle(x=1.0, y=1.0; periodic = false)
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).
"""
function billiard_rectangle(x=1.0, y=1.0; periodic = false)

  bt = Obstacle[]
  if periodic == false
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
    leftw.partner = rightw
    rightw.partner = leftw
    sp = [o,y]; ep = [x, y]; n = [o,-y]
    topw = PeriodicWall(sp, ep, n, "Top periodic boundary")
    sp = [o,o]; ep = [x, o]; n = [o,y]
    botw = PeriodicWall(sp, ep, n, "Bottom periodic boundary")
    topw.partner = botw
    botw.partner = topw
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
  bt = billiard_rectangle(x,y; periodic = periodic)
  if periodic && r>=x/2 || r>=y/2
    es = "Disk radius too big for a periodic Sinai billiard.\n"
    es*= "Obstacles cannot overlap with `PeriodicWall`s."
    error(es)
  end
  c = [x/2, y/2]
  centerdisk = Disk(c, r, "Disk")
  push!(bt, centerdisk)
end
