using StaticArrays
using PyPlot
import Base.show

# Find a way to find the next collision in an "infinite lattice", just like DP Sanders
# did for straight propagation. Do not know if it is possible.

####################################################
## Particles
####################################################
"""
Particle supertype.
"""
abstract AbstractParticle{T<:AbstractFloat}


"""
    Particle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table.
# Fields:
* `pos::SVector{2,T}` : Current position vector, always within the the limits of `cellsize`.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at.
(Used only in periodic billiards)
"""
type Particle{T<:AbstractFloat} <: AbstractParticle{T}
  pos::SVector{2,T}
  vel::SVector{2,T}
  current_cell::SVector{2,T}
end
"""
    Particle{T}(ic::Vector{T})
Constructor accepting initial conditions `[x0, y0, φ0]`.
"""
function Particle{T<:AbstractFloat}(ic::Vector{T})
  φ0 = ic[3]
  pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
  Particle(pos, vel, SVector{2,T}(zero(T), zero(T)))
end
Particle{T<:AbstractFloat}(x0::T, y0::T, φ0::T) = Particle([x0, y0, φ0])
Particle() = Particle(rand(), rand(), rand()*2π)


"""
    MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table with perpendicular magnetic field.
# Fields:
* `pos::SVector{2,T}` : Current position vector, always within the the limits of `cellsize`.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at.
(Used only in periodic billiards)
* `ω::T` : Angular velocity of rotational motion. Radius of rotation is `r=1/ω`
"""
type MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
  pos::SVector{2,T}
  vel::SVector{2,T}
  current_cell::SVector{2,T}
  ω::T
  function MagneticParticle(pos::SVector{2,T}, vel::SVector{2,T},
                            current_cell::SVector{2,T}, ω::T)
    if ω==0
      error("Angular frequency cannot be 0.")
    end
    new(pos, vel, current_cell, ω)
  end
end
"""
    MagneticParticle{T}(ic::Vector{T}, ω)
Constructor accepting initial conditions `[x0, y0, φ0]`.
"""
function MagneticParticle{T<:AbstractFloat}(ic::Vector{T}, ω::T)
  φ0 = ic[3]
  pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
  MagneticParticle(pos, vel, SVector{2,T}(zero(T), zero(T)), ω)
end
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)

function Particle{T<:AbstractFloat}(p::MagneticParticle{T}; use_cell = true)
  pos = p.pos
  cell = ifelse(use_cell, p.current_cell, SVector{T,2}(0,0))
  Particle(pos, p.vel, cell)
end




"""
    cyclotron{T<:AbstractFloat}(p::AbstractParticle{T}, ω::T)
Return center and radius of circular motion performed by the particle based on
`p.pos` and `p.vel`.
"""
function cyclotron{T<:AbstractFloat}(ω::T, p::AbstractParticle{T})
  c::SVector{2, T} = p.pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  return c, r
end

"""
    cyclotron{T}(p::MagneticParticle{T}, use_cell = false)
Return center and radius of circular motion performed by the particle based on
`p.pos` (or `p.pos + p.current_cell`) and `p.vel`.
"""
function cyclotron{T<:AbstractFloat}(p::MagneticParticle{T}, use_cell = false)
  ω = p.ω
  pos = ifelse(use_cell, p.pos + p.current_cell, p.pos)
  c::SVector{2, T} = pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  return c, r
end

####################################################
## Obstacles
####################################################
"""
Abstract Type, only to serve as a node in Type graph.
"""
abstract Obstacle{T<:AbstractFloat}
abstract Circular{T<:AbstractFloat} <: Obstacle{T}
immutable Disk{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    name::String
    #Inner constructor, do not add {T} after name!
    function Disk(c::SVector{2,T}, r::T, name::String)
      new(c, abs(r), name)
    end
end
function Disk{T<:AbstractFloat}(center::Vector{T}, radius::T, name::String)
  Disk{T}(SVector{2, T}(center), radius, name)
end
immutable Circle{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    name::String
    #Inner constructor, do not add {T} after name!
    function Circle(c::SVector{2,T}, r::T, name::String)
      new(c, abs(r), name)
    end
end
function Circle{T<:AbstractFloat}(center::Vector{T}, radius::T, name::String)
  Circle{T}(SVector{2, T}(center), radius, name)
end

type Antidot{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    inside::Bool
    name::String
    #Inner constructor, do not add {T} after name!
    function Antidot(c::SVector{2,T}, r::T, inside::Bool, name::String)
      new(c, abs(r), inside, name)
    end
end
function Antidot{T<:AbstractFloat}(center::Vector{T}, radius::T, inside::Bool, name::String)
  Antidot{T}(SVector{2, T}(center), radius, inside, name)
end
function Antidot{T<:AbstractFloat}(center::Vector{T}, radius::T, name::String)
  Antidot{T}(SVector{2, T}(center), radius, false, name)
end

abstract Wall{T<:AbstractFloat} <: Obstacle{T}
"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
A wall imposing specular reflection during collision.
# Fields
* `sp::SVector{2,T}` : Starting point of the Wall (non-negative).
* `ep::SVector{2,T}` : Ending point of the Wall (non-negative).
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the particle *will
come from* (to the inside the billiard table). The size of the vector is irrelevant, but
must be normal to the wall.
* `name::String` : The name of the object, e.g. "left wall".
"""
immutable FiniteWall{T<:AbstractFloat} <: Wall{T}
  sp::SVector{2,T}
  ep::SVector{2,T}
  normal::SVector{2,T}
  name::String
  #Inner constructor, do not add {T} after name
  function FiniteWall(sp::SVector{2,T}, ep::SVector{2,T},
                      normal::SVector{2,T}, name::String)
    d = dot(normal, ep-sp)
    if d != zero(T)
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, normal, name)
  end
end
function FiniteWall{T<:AbstractFloat}(sp::Vector{T}, ep::Vector{T},
  n::Vector{T}, name::String)
  # Be sure you have {T} after the name of the Constructor call because the type is called.
  FiniteWall{T}(SVector{2, T}(sp), SVector{2, T}(ep), SVector{2, T}(n), name)
end
function NullWall{T<:AbstractFloat}(::Type{T})
  a = 1000.0*one(T)
  b = one(T)
  o = zero(T)
  FiniteWall([a, o], [a, b],  [-b, o], "Null wall")
end

#pretty print
show{T}(io::IO, w::FiniteWall{T}) =
    print(io, "$(w.name) {$T}\n",
    "start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")

"""
    PeriodicWall{T<:AbstractFloat} <: Wall{T}
This type of Wall imposes periodic boundary conditions during collision.
# Fields
* `sp::SVector{2,T}` : Starting point of the Wall (non-negative).
* `ep::SVector{2,T}` : Ending point of the Wall (non-negative).
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the particle *will
come from* (to the inside the billiard table). The size of the vector is **important**.
This vector is added to a particle's `pos` during collision. Therefore the size of the
normal vector must be correctly associated with the size of the cell (field `cellsize`) of
the `PeriodicParticle` instance.
* `partner::PeriodicWall{T}` : Associated periodic partner of a PeriodicWall. This field is
undefined during instantiation. It must be defined during the construction of the Billiard
Table.
"""
type PeriodicWall{T<:AbstractFloat} <: Wall{T}
  sp::SVector{2,T}
  ep::SVector{2,T}
  normal::SVector{2,T}
  name::String
  partner::PeriodicWall{T}
  function PeriodicWall(sp, ep, normal, name)
    for el in vcat(sp, ep)
      el < 0 && error("Negative values are not allowed for periodic wall objects.")
    end
    d = abs(dot(normal, ep-sp))
    if d != zero(T)
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, normal, name) # The `partner` field is uninitialized.
  end
end
"""
    PeriodicWall{T<:AbstractFloat}(sp::Vector{T}, ep::Vector{T}, n::Vector{T})
Constructor accepting 3 `Vector{T}` arguments.
"""
function PeriodicWall{T<:AbstractFloat}(sp::Vector{T}, ep::Vector{T},
  n::Vector{T}, name::String)
  # I have to call the TYPE here. That is why I NEED TO PUT {T}! For functions, I MUST NOT
  # put {T}!
  PeriodicWall{T}(SVector{2, T}(sp), SVector{2, T}(ep), SVector{2, T}(n), name)
end
show{T}(io::IO, w::PeriodicWall{T}) =
    print(io, "$(w.name) {$T}\n",
    "start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)\n",
    "periodic partner: $(w.partner.name)")

#normalvec will be only called internally
normalvec{T<:AbstractFloat}(disk::Circular{T}, pos::SVector{2,T}) = normalize(pos - disk.c)
normalvec{T<:AbstractFloat}(wall::Wall{T}, pos::SVector{2,T}) = normalize(wall.normal)

####################################################
## Resolve Collisions
####################################################
#All these functions will become dependend on velocity angle for ray-splitting billiards.
function specular!{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T})
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
end

function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::Obstacle{T},
  printstuff = false)

  #println("Resolving collision with $(o.name)")
  #id = randstring(4)

  dist = distance(p, o)
  if abs(dist) > 1e-4
    println("In resolve with $(o.name), dist = $dist")
    error("Distance unreasonably big")
  end
  dt = 0.0
  #printstuff = false
  if dist < 0.0
    if printstuff
      println("In resolve with $(o.name), dist = $dist")
    end
    #dt = 10lct(p, o, dist) #linearized collision time, should be negative
    dt = 10lct(p, o, dist) #linearized collision time, should be negative
    printstuff && println("2*l.c.t. = $dt")
    if dt > 0
      n = normalvec(o, p.pos)
      println("dot(p.vel, n) = $(dot(p.vel, n))")
      println("(should be negative)")
      error("lct should be negative")
    end
    # Check the orders of magnitude:
    vxt=p.vel[1]*dt; vyt= p.vel[2]*dt
    #exponent of velocity*timeL
    ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
    #exponent of position
    pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
    # Ensure that the position will be changed
    if pe - ve > 16.0
      dt *= 10^(pe - ve - 16.0)
    end
    #println("dt = $dt\n")
    # Propagate backwards:
    p.pos += [p.vel[1]*dt, p.vel[2]*dt]
  elseif printstuff
    println("Good distance for $(o.name), dist = $dist")
  end

  #TEST TO BE SURE
  if distance(p, o)<0.0
    error("Dist = $(distance(p,o)) after back-propagation")
  end


  # Perform specular reflection:
  n = normalvec(o, p.pos)
  p.vel = p.vel - 2*dot(n, p.vel)*n
  printstuff && println("-----------------------------------------------------\n")
  return dt
end

function resolvecollision!{T<:AbstractFloat}(p::AbstractParticle{T}, o::PeriodicWall{T},
  printstuff = false)

    #println("Resolving collision with $(o.name)")
    #id = randstring(4)
    dist = distance(p, o)
    if abs(dist) > 1e-4
      println("In resolve with $(o.name), dist = $dist")
      error("Distance unreasonably big")
    end
    t = 0.0
    #printstuff = false

    if dist > 0.0
      if printstuff

        println("In resolve with $(o.name), dist = $dist")
      end
      t = 10lct(p, o, dist) #linearized collision time, should be positive
      if t < 0
        n = normalvec(o, p.pos)
        println("dot(p.vel, n) = $(dot(p.vel, n))")
        println("(should be negative)")
        error("lct should be positive")
      end
      # Check the orders of magnitude:
      vxt=p.vel[1]*t; vyt= p.vel[2]*t
      #exponent of velocity*timeL
      ve = exponent(min(abs(vxt), abs(vyt)))*0.3010299956639812
      #exponent of position
      pe = exponent(max(abs(p.pos[1]), abs(p.pos[2])))*0.3010299956639812
      # Ensure that the position will be changed
      if pe - ve > 16
        t *= 10^(pe - ve - 16)
      end
      # Propagate forwards:
      p.pos += [p.vel[1]*t, p.vel[2]*t]
    elseif printstuff
      println("Good distance for $(o.name), dist = $dist")
    end
  #TEST TO BE SURE
  if distance(p, o)>0.0
    error("Dist = $(distance(p,o)) after forward-propagation")
  end

  #perform periodicity
  p.pos += o.normal
  p.current_cell -= o.normal
  printstuff && println("-----------------------------------------------------\n")
  return t
end

####################################################
## Distances
####################################################
"""
    distance(p::AbstractParticle{T}, o::Obstacle{T})
Return the **signed** distance between Particle and Obstacle, based on `p.pos`.
Positive distance corresponds to the particle being inside the *allowed* region
of the Billiard Table/Obstacle.
    distance(p::AbstractParticle{T}, bt::Vector{Obstacle{T}})
Return minimum `distance(p, obst)` for all `obst` in `bt`.
"""
function distance{T<:AbstractFloat}(p::AbstractParticle{T}, bt::Vector{Obstacle{T}})
  mindist = convert(T, Inf)
  for obst in bt
    dist = distance(p, obst)
    if dist <= mindist; mindist = dist; end
  end
  return mindist
end

function distance{T<:AbstractFloat}(p::AbstractParticle{T}, w::Wall{T})
  v1 = p.pos - w.sp
  dot(v1, w.normal)
end

function distance{T<:AbstractFloat}(p::AbstractParticle{T}, d::Disk{T})
  norm(p.pos - d.c) - d.r
end

function distance{T<:AbstractFloat}(p::AbstractParticle{T}, d::Circle{T})
  d.r - norm(p.pos - d.c)
end

####################################################
## Initial Conditions
####################################################
"""
    cellsize{T<:AbstractFloat}(bt::Vector{Obstacle{T}})
Return the delimiters `xmin, ymin, xmax, ymax` of the "cell" that is defined by
the given billiard table. Used in `randominside()` and error checking.
"""
function cellsize{T<:AbstractFloat}(bt::Vector{Obstacle{T}})

  xmin = ymin = T(Inf)
  xmax = ymax = T(-Inf)
  #test if there is Circle in bt:
  if any(x -> isa(x, Circle), bt)
    i = find(x -> isa(x, Circle), bt)
    for j in i
      xmin2 = bt[j].c[1] -  bt[j].r
      ymin2 = bt[j].c[2] -  bt[j].r
      xmax2 = bt[j].c[1] +  bt[j].r
      ymax2 = bt[j].c[2] +  bt[j].r

      xmin = min(xmin, xmin2)
      ymin = min(ymin, ymin2)
      xmax = max(xmax, xmax2)
      ymax = max(ymax, ymax2)
    end
    return xmin, ymin, xmax, ymax
  end
  #Else, use the walls:
  for obst in bt
    if typeof(obst) <: Wall
      xmin2 = min(obst.sp[1], obst.ep[1])
      ymin2 = min(obst.sp[2], obst.ep[2])
      xmax2 = max(obst.sp[1], obst.ep[1])
      ymax2 = max(obst.sp[2], obst.ep[2])
      xmin = min(xmin, ymin2)
      ymin = min(ymin, ymin2)
      xmax = max(xmax, xmax2)
      ymax = max(ymax, ymax2)
    end
    if any(x -> isa(x, PeriodicWall), bt) && [xmin, ymin] != [0,0]
      error("Periodic billiard tables must have (0,0) as the bottom left corner.")
    end
  end
  return xmin, ymin, xmax, ymax
end

function randominside{T<:AbstractFloat}(bt::Vector{Obstacle{T}})

  xmin, ymin, xmax, ymax = cellsize(bt)
  f = convert(T, rand())
  while f == 0 || f==1/4 || f==1/2 || f == 3/4
    f = convert(T, rand())
  end
  φ0 = f*2π

  xp = rand()*(xmax-xmin) + xmin
  yp = rand()*(ymax-ymin) + ymin
  p = Particle(xp, yp, φ0)

  dist = distance(p, bt)
  while dist <= 0.0

    xp = rand()*(xmax-xmin) + xmin
    yp = rand()*(ymax-ymin) + ymin
    p.pos = [xp, yp]
    dist = distance(p, bt)
  end

  p = Particle(xp, yp, φ0)
end




####################################################
## Famous/Standard Billiards
####################################################
function billiard_rectangle{T<:AbstractFloat}(x::T, y::T)

  bt = Obstacle{T}[]
  o = zero(T)
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
billiard_rectangle() = billiard_rectangle(1.0,1.0)

function billiard_sinai{T<:AbstractFloat}(r::T, x::T=one(T), y::T=one(T))
  bt = billiard_rectangle(x,y)
  c = [x/2, y/2]
  centerdisk = Disk(c, r, "Disk")
  push!(bt, centerdisk)
end

function billiard_sinai_periodic{T<:AbstractFloat}(r::T, x::T=one(T), y::T=one(T))
  bt = Obstacle{T}[]
  o = zero(T)

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
