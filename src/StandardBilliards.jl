####################################################
## Famous/Standard Billiards
####################################################
"""
    billiard_rectangle(x=1.0, y=1.0)
Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).
"""
function billiard_rectangle(x=1.0, y=1.0)

  bt = Obstacle[]
  o = 0.0
  sp = [o,o]; ep = [o, y]; n = [x,o]
  leftw = FiniteWall(sp, ep, n, "Left wall")
  sp = [x,o]; ep = [x, y]; n = [-x,o]
  rightw = FiniteWall(sp, ep, n, "Right wall")
  sp = [o,y]; ep = [x, y]; n = [o,-y]
  topw = FiniteWall(sp, ep, n, "Top wall")
  sp = [o,o]; ep = [x, o]; n = [o,y]
  botw = FiniteWall(sp, ep, n, "Bottom wall")
  push!(bt, leftw, rightw, topw, botw)
end

"""
    billiard_sinai(r, x=1.0, y=1.0)
Return a vector of obstacles that defines a closed Sinai billiard of size (`x`, `y`) with
a disk in its center, of radius `r`.
"""
function billiard_sinai(r, x=1.0, y=1.0)
  bt = billiard_rectangle(x,y)
  c = [x/2, y/2]
  centerdisk = Disk(c, r, "Disk")
  push!(bt, centerdisk)
end

"""
    billiard_sinai_periodic(r, x=1.0, y=1.0)
Return a vector of obstacles that defines a periodic Sinai billiard (aka Lorentz Gas)
of size (`x`, `y`) with a disk in its center, of radius `r`.
"""
function billiard_sinai_periodic(r, x=1.0, y=1.0)
  bt = Obstacle[]
  o = 0.0

  if r>=x/2 || r>=y/2
    error("Disk radius too big for a periodic Sinai billiard.")
  end

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
  c = [x/2, y/2]
  centerdisk = Disk(c, r, "Disk")
  push!(bt, centerdisk)
end
