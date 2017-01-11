using StaticArrays
import Base.show

####################################################
## Particles
####################################################
"""
    AbstractParticle{T<:AbstractFloat}
Particle supertype.
"""
abstract AbstractParticle{T<:AbstractFloat}

"""
    Particle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table.
# Fields:
* `pos::SVector{2,T}` : Current position vector.
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
* `pos::SVector{2,T}` : Current position vector.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at.
(Used only in periodic billiards)
* `omega::T` : Angular velocity (cyclic frequency) of rotational motion.
Radius of rotation is `r=1/omega`.
"""
type MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
  pos::SVector{2,T}
  vel::SVector{2,T}
  current_cell::SVector{2,T}
  omega::T
  function MagneticParticle(pos::SVector{2,T}, vel::SVector{2,T},
                            current_cell::SVector{2,T}, omega::T)
    if omega==0
      error("Angular frequency cannot be 0.")
    end
    new(pos, vel, current_cell, omega)
  end
end
"""
    MagneticParticle{T<:AbstractFloat}(ic::Vector{T}, omega::T)
Constructor accepting initial conditions `[x0, y0, φ0]`.
"""
function MagneticParticle{T<:AbstractFloat}(ic::Vector{T}, ω::T)
  φ0 = ic[3]
  pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
  return MagneticParticle{T}(pos, vel, SVector{2,T}(zero(T), zero(T)), ω)
end
MagneticParticle{T<:AbstractFloat}(x0::T, y0::T, φ0::T, ω::T) =
MagneticParticle([x0, y0, φ0], ω)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)

"""
    magnetic2standard(p::MagneticParticle; use_cell = true)
Create a standard particle from a `MagneticParticle`.
"""
function magnetic2standard{T<:AbstractFloat}(p::MagneticParticle{T}; use_cell = true)
  pos = p.pos
  cell = ifelse(use_cell, p.current_cell, SVector{T,2}(0,0))
  Particle(pos, p.vel, cell)
end

"""
    magnetic2standard(p::Particle, omega; use_cell = true)
Create a magnetic particle from a `Particle`.
"""
function standard2magnetic{T<:AbstractFloat}(p::Particle{T}, ω::T; use_cell = true)
  pos = p.pos
  cell = ifelse(use_cell, p.current_cell, SVector{T,2}(0,0))
  MagneticParticle(pos, p.vel, cell, ω)
end


"""
    cyclotron(omega, p::AbstractParticle)
Return center and radius of circular motion performed by the particle based on
`p.pos` and `p.vel`.
"""
function cyclotron{T<:AbstractFloat}(ω::T, p::AbstractParticle{T})
  c::SVector{2, T} = p.pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  return c, r
end

"""
    cyclotron(p::MagneticParticle, use_cell = false)
Return center and radius of circular motion performed by the particle based on
`p.pos` (or `p.pos + p.current_cell`) and `p.vel`.
"""
function cyclotron{T<:AbstractFloat}(p::MagneticParticle{T}, use_cell = false)
  ω = p.omega
  pos = ifelse(use_cell, p.pos + p.current_cell, p.pos)
  c::SVector{2, T} = pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  return c, r
end

####################################################
## Obstacles
####################################################
"""
    Obstacle{T<:AbstractFloat}
Obstacle supertype.
"""
abstract Obstacle{T<:AbstractFloat}
"""
    Circular{T<:AbstractFloat} <: Obstacle{T}
Circular obstacle supertype.
"""
abstract Circular{T<:AbstractFloat} <: Obstacle{T}

"""
    Disk{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle with propagation allowed outside of the circle.
# Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `name::String` : Some name given for user convenience.
"""
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

"""
    Circle{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle with propagation allowed inside of the circle.
# Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `name::String` : Some name given for user convenience.
"""
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

"""
    Antidot{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that allows propagation both inside and outside of the Circle.
Used in ray-splitting billiards.
# Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `inside::Bool` : Flag that keeps track whether the particle is inside or outside
of the disk.
* `name::String` : Name of the obstacle given for user convenience.
"""
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

"""
    Wall{T<:AbstractFloat} <: Obstacle{T}
Wall obstacle supertype.
"""
abstract Wall{T<:AbstractFloat} <: Obstacle{T}

"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision.
# Fields
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the particle *will
come from before a collision* (pointing towards the inside the billiard table).
The size of the vector is irrelevant.
* `name::String` : Name of the obstacle, e.g. "left wall", given for user convenience.
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
Wall obstacle that imposes periodic boundary conditions during collision.
# Fields
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the particle *will
come from* (to the inside the billiard table). The size of the vector is **important**.
This vector is added to a particle's `pos` during collision. Therefore the size of the
normal vector must be correctly associated with the size of the cell. Also, the partners
(see field `partner`) must have same size in their normal vectors.
* `name::String` : Name of the obstacle, e.g. "left boundary", given for user convenience.
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
"""
    normalvec(obst::Obstacle{T}, position::SVector{2,T})
Return the vector normal to the obstacle from the current particle position (which is
assumed to be on top of the obstacle's boundary).
"""
normalvec{T<:AbstractFloat}(disk::Circular{T}, pos::SVector{2,T}) = normalize(pos - disk.c)
normalvec{T<:AbstractFloat}(wall::Wall{T}, pos::SVector{2,T}) = normalize(wall.normal)


####################################################
## Distances
####################################################
# This should not be exported
"""
    distance(p::AbstractParticle, o::Obstacle)
Return the **signed** distance between particle `p` and obstacle `o`, based on `p.pos`.
Positive distance corresponds to the particle being inside the *allowed* region
of the Obstacle.
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


"""
    randominside(bt::Vector{Obstacle})
    randominside(bt::Vector{Obstacle}, omega)
Return a particle with correct (allowed) initial conditions inside the given billiard
table defined by the vector `bt`. If supplied with a second argument the type of
the returned particle is `MagneticParticle`, with angular velocity `omega`.
Else, it is `Particle`.
"""
function randominside{T<:AbstractFloat}(bt::Vector{Obstacle{T}})

  xmin, ymin, xmax, ymax = cellsize(bt)
  f = convert(T, rand())
  while f == 0 || f==1/4 || f==1/2 || f == 3/4
    f = convert(T, rand())
  end
  φ0 = f*2π

  xp = rand()*(xmax-xmin) + xmin
  yp = rand()*(ymax-ymin) + ymin
  p = Particle([xp, yp, φ0])

  dist = distance(p, bt)
  while dist <= 0.0

    xp = rand()*(xmax-xmin) + xmin
    yp = rand()*(ymax-ymin) + ymin
    p.pos = [xp, yp]
    dist = distance(p, bt)
  end

  return p
end

function randominside{T<:AbstractFloat}(ω::T, bt::Vector{Obstacle{T}})

  xmin, ymin, xmax, ymax = cellsize(bt)
  f = convert(T, rand())
  while f == 0 || f==1/4 || f==1/2 || f == 3/4
    f = convert(T, rand())
  end
  φ0 = f*2π

  xp = rand()*(xmax-xmin) + xmin
  yp = rand()*(ymax-ymin) + ymin
  p = MagneticParticle([xp, yp, φ0], ω)

  dist = distance(p, bt)
  while dist <= 0.0

    xp = rand()*(xmax-xmin) + xmin
    yp = rand()*(ymax-ymin) + ymin
    p.pos = [xp, yp]
    dist = distance(p, bt)
  end

  return p
end
