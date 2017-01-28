using StaticArrays
import Base.show
export AbstractParticle, Particle, MagneticParticle, magnetic2standard,
standard2magnetic, cyclotron, Obstacle, Disk, Antidot, FiniteWall, PeriodicWall,
normalvec, randominside, SplitterWall, distance, cellsize

####################################################
## Particles
####################################################
"""
    AbstractParticle
Particle supertype.
"""
abstract AbstractParticle

"""
    Particle <: AbstractParticle
Two-dimensional particle in a billiard table.
## Fields:
* `pos::SVector{2,Float64}` : Current position vector.
* `vel::SVector{2,Float64}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,Float64}` : Current "cell" the particle is located at.
  (Used only in periodic billiards)
## Additional constructors:
```julia
Particle{T<:Real}(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x::Real, y::Real, φ::Real)
Particle() = Particle(rand(), rand(), rand()*2π)
```
"""
type Particle <: AbstractParticle
  pos::SVector{2,Float64}
  vel::SVector{2,Float64}
  current_cell::SVector{2,Float64}
end

function Particle{T<:Real}(icr::Vector{T})
  ic = Float64.(icr)
  φ0 = ic[3]
  pos = SVector{2,Float64}(ic[1:2]); vel = SVector{2,Float64}(cos(φ0), sin(φ0))
  return Particle(pos, vel, SVector{2,Float64}(0,0))
end
Particle(x0::Real, y0::Real, φ0::Real) = Particle([x0, y0, φ0])
Particle() = Particle(rand(), rand(), rand()*2π)
show(io::IO, p::Particle) =
    print(io, "MagneticParticle\n",
    "position: $(p.pos)\nvelocity: $(p.vel)")

"""
    MagneticParticle <: AbstractParticle
Two-dimensional particle in a billiard table with perpendicular magnetic field.
## Fields:
* `pos::SVector{2,Float64}` : Current position vector.
* `vel::SVector{2,Float64}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,Float64}` : Current "cell" the particle is located at.
  (Used only in periodic billiards)
* `omega::Float64` : Angular velocity (cyclic frequency) of rotational motion.
  Radius of rotation is `r=1/omega`.
## Additional constructors:
```julia
MagneticParticle{T<:Real}(ic::Vector{T}, ω::Real) #where ic = [x0, y0, φ0]
MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)
```
"""
type MagneticParticle <: AbstractParticle
  pos::SVector{2,Float64}
  vel::SVector{2,Float64}
  current_cell::SVector{2,Float64}
  omega::Float64
  function MagneticParticle(pos::SVector{2,Float64}, vel::SVector{2,Float64},
                            current_cell::SVector{2,Float64}, omega::Float64)
    if omega==0
      error("Angular velocity cannot be 0.")
    end
    new(pos, vel, current_cell, float(omega))
  end
end

function MagneticParticle{T<:Real}(ic::Vector{T}, ω::Real)
  φ0 = ic[3]
  pos = SVector{2,Float64}(ic[1:2]); vel = SVector{2,Float64}(cos(φ0), sin(φ0))
  return MagneticParticle(pos, vel, SVector{2,Float64}(0,0), float(ω))
end
MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real) =
MagneticParticle([x0, y0, φ0], ω)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)

show(io::IO, p::MagneticParticle) =
    print(io, "MagneticParticle\n",
    "position: $(p.pos)\nvelocity: $(p.vel)\nω: $(p.omega)")

"""
```julia
magnetic2standard(p::MagneticParticle, use_cell = true)
```
Create a standard `Particle` from a `MagneticParticle`.
"""
function magnetic2standard(p::MagneticParticle, use_cell = true)
  pos = p.pos
  cell = use_cell ? current_cell : SVector{Float64,2}(0,0)
  return Particle(pos, p.vel, cell)
end

"""
```julia
standard2magnetic(omega, p::Particle, use_cell = true)
```
Create a `MagneticParticle` from a `Particle`.
"""
function standard2magnetic(ω::Real, p::Particle, use_cell = true)
  pos = p.pos
  cell = use_cell ? current_cell : SVector{Float64,2}(0,0)
  return MagneticParticle(pos, p.vel, cell, ω)
end


"""
```julia
cyclotron(p::MagneticParticle, use_cell = false)
```
Return center and radius of circular motion performed by the particle based on
`p.pos` (or `p.pos + p.current_cell`) and `p.vel`.
"""
function cyclotron(p::MagneticParticle, use_cell = false)
  ω = p.omega
  pos = ifelse(use_cell, p.pos + p.current_cell, p.pos)
  c::SVector{2, Float64} = pos - (1/ω)*[p.vel[2], -p.vel[1]]
  r = abs(1/ω)
  return c, r
end

####################################################
## Obstacles
####################################################
"""
    Obstacle
Obstacle supertype.
"""
abstract Obstacle
"""
    Circular <: Obstacle
Circular obstacle supertype.
"""
abstract Circular <: Obstacle

"""
    Disk <: Circular
Disk-like obstacle with propagation allowed outside of the circle.
## Fields:
* `c::SVector{2,Float64}` : Center.
* `r::Float64` : Radius.
* `name::String` : Some name given for user convenience.
Constructors accept any vectors convertible to SVector{2,Float64}.
"""
immutable Disk <: Circular
    c::SVector{2,Float64}
    r::Float64
    name::String
    function Disk(c, r::Real, name::String = "Disk")
      new(SVector{2,Float64}(c), abs(r), name)
    end
end

"""
    Antidot <: Circular
Disk-like obstacle that allows propagation both inside and outside of the disk.
Used in ray-splitting billiards.
## Fields:
* `c::SVector{2,Float64}` : Center.
* `r::Float64` : Radius.
* `where::Bool` : Flag that keeps track of where the particle is currently
  propagating. `true` stands for *outside* the disk, `false` for *inside* the disk.
* `name::String` : Name of the obstacle given for user convenience.
Constructors accept any vectors convertible to SVector{2,Float64}.
"""
type Antidot <: Circular
    c::SVector{2,Float64}
    r::Float64
    where::Bool
    name::String
    function Antidot(c, r::Real, where::Bool = true, name::String = "Antidot")
      new(SVector{2,Float64}(c), abs(r), where, name)
    end
end
function Antidot(c, r::Real, name::String)
  Antidot(SVector{2,Float64}(c), abs(r), true, name)
end

show(io::IO, w::Circular) =
    print(io, "$(w.name)\n",
    "center: $(w.c)\nradius: $(w.r)")


"""
    Wall <: Obstacle
Wall obstacle supertype.
"""
abstract Wall <: Obstacle

"""
    FiniteWall <: Wall
Wall obstacle imposing specular reflection during collision.
## Fields:
* `sp::SVector{2,Float64}` : Starting point of the Wall.
* `ep::SVector{2,Float64}` : Ending point of the Wall.
* `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside the
  billiard table). The size of the vector is irrelevant.
* `name::String` : Name of the obstacle, e.g. "left wall", given for user convenience.
Constructors accept any vectors convertible to SVector{2,Float64}.
"""
immutable FiniteWall <: Wall
  sp::SVector{2,Float64}
  ep::SVector{2,Float64}
  normal::SVector{2,Float64}
  name::String
  #Inner constructor, do not add {Float64} after name
  function FiniteWall(sp::SVector{2,Float64}, ep::SVector{2,Float64},
    normal::SVector{2,Float64}, name::String = "Wall")

    n = normalize(normal)
    d = dot(n, ep-sp)
    if abs(d) > 1e-15
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, n, name)
  end
end
function FiniteWall(sp, ep, n, name::String = "wall")
  FiniteWall(SVector{2, Float64}(sp), SVector{2, Float64}(ep),
  SVector{2, Float64}(n), name)
end

#pretty print
show(io::IO, w::FiniteWall) =
    print(io, "$(w.name)\n",
    "start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")

"""
    PeriodicWall <: Wall
Wall obstacle that imposes periodic boundary conditions upon collision.
## Fields:
* `sp::SVector{2,Float64}` : Starting point of the Wall.
* `ep::SVector{2,Float64}` : Ending point of the Wall.
* `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the
  particle *will come from* (to the inside the billiard table).
  The size of the vector is **important**.
  This vector is added to a particle's `pos` during collision. Therefore the
  size of the normal vector must be correctly associated with the size of the
  periodic cell.
* `name::String` : Name of the obstacle, e.g. "left boundary", given for user
  convenience.
Constructors accept any vectors convertible to SVector{2,Float64}.
"""
immutable PeriodicWall <: Wall
  sp::SVector{2,Float64}
  ep::SVector{2,Float64}
  normal::SVector{2,Float64}
  name::String
  function PeriodicWall(sp::SVector{2,Float64}, ep::SVector{2,Float64},
                        normal::SVector{2,Float64}, name::String)
    d = dot(normal, ep-sp)
    if abs(d) > 1e-15
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, normal, name)
  end
end

function PeriodicWall(sp, ep, n, name::String = "Periodic wall")
  return PeriodicWall(SVector{2, Float64}(sp), SVector{2, Float64}(ep),
  SVector{2, Float64}(n), name)
end

show(io::IO, w::PeriodicWall) =
print(io, "$(w.name)\n",
"start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")

"""
    SplitterWall <: Wall
Wall obstacle imposing specular reflection during collision.
## Fields:
* `sp::SVector{2,Float64}` : Starting point of the Wall.
* `ep::SVector{2,Float64}` : Ending point of the Wall.
* `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision*.
  The size of the vector is irrelevant.
* `where::Bool` : Flag that keeps track of where the particle is currently
  propagating. `true` is associated with the `normal` vector the wall is
  instantiated with.
* `name::String` : Name of the obstacle, e.g. "ray-splitting wall 1",
  given for user convenience.
Constructors accept any vectors convertible to SVector{2,Float64}.
"""
type SplitterWall <: Wall
  sp::SVector{2,Float64}
  ep::SVector{2,Float64}
  normal::SVector{2,Float64}
  where::Bool
  name::String
  #Inner constructor, do not add {Float64} after name
  function SplitterWall(sp::SVector{2,Float64}, ep::SVector{2,Float64},
    normal::SVector{2,Float64}, where::Bool = true, name::String = "Wall")

    n = normalize(normal)
    d = dot(n, ep-sp)
    if abs(d) > 1e-15
      error("Normal vector is not actually normal to the wall")
    end
    new(sp, ep, n, where, name)
  end
end

function SplitterWall(sp, ep, n, where::Bool=true, name::String="Ray-splitting wall")
  return SplitterWall(SVector{2, Float64}(sp), SVector{2, Float64}(ep),
  SVector{2, Float64}(n), where, name)
end

show(io::IO, w::SplitterWall) =
print(io, "$(w.name)\n",
"start point: $(w.sp)\nend point: $(w.ep)\n",
"current normal: $((2*Int(w.where)- 1)*w.normal)")

"""
```julia
normalvec(obst::Obstacle, position)
```
Return the vector normal to the obstacle at the given position (which is
assumed to be very close to the obstacle's boundary).
The normal vector of any Obstacle must be looking towards
the direction a particle is expected to come from.
"""
normalvec(wall::FiniteWall, pos) = wall.normal
normalvec(w::PeriodicWall, pos) = normalize(w.normal)
normalvec(w::SplitterWall, pos) = (2*Int(w.where)- 1)*w.normal
normalvec(disk::Circular, pos) = normalize(pos - disk.c)
normalvec(a::Antidot, pos) = (2*Int(a.where)- 1)*normalize(pos - a.c)

####################################################
## Distances
####################################################
"""
```julia
distance(p::AbstractParticle, o::Obstacle)
```
Return the **signed** distance between particle `p` and obstacle `o`, based on
`p.pos`. Positive distance corresponds to the particle being on the *allowed* region
of the Obstacle. E.g. for a `Disk`, the distance is positive when the particle is
outside of the disk, negative otherwise.
```julia
distance(p::AbstractParticle, bt::Vector{Obstacle})
```
Return minimum `distance(p, obst)` for all `obst` in `bt`, which can be negative.
"""
function distance(p::AbstractParticle, bt::Vector{Obstacle})
  mindist = Inf
  for obst in bt
    dist = distance(p, obst)
    if dist <= mindist; mindist = dist; end
  end
  return mindist
end

# function distance(p::AbstractParticle, w::Wall)
#   v1 = p.pos - w.sp
#   dot(v1, w.normal)
# end

function distance(p::AbstractParticle, w::Wall)
  v1 = p.pos - w.sp
  dot(v1, normalvec(w, p.pos))
end

# no new distance needed for SplitterWall because the `where` field
# has the necessary information to give the correct dinstance.
# function distance(p::AbstractParticle, w::SplitterWall)
#   v1 = p.pos - w.sp
#   (2*Int(a.where)- 1)*dot(v1, w.normal)
# end

distance(p::AbstractParticle, d::Disk) = norm(p.pos - d.c) - d.r

function distance(p::AbstractParticle, a::Antidot)
  (2*Int(a.where)- 1)*(norm(p.pos - a.c) - a.r)
end


####################################################
## Initial Conditions
####################################################

"""
```julia
    cellsize(obst::Obstacle)
    cellsize(bt::Vector{Obstacle})
```
Return the delimiters `xmin, ymin, xmax, ymax` of the given obstacle/billiard table.

Used in `randominside()` and error checking.
"""
function cellsize(d::Disk)
  xmin = ymin = Inf
  xmax = ymax = -Inf
  return xmin, ymin, xmax, ymax
end

function cellsize(w::Wall)
  xmin = min(w.sp[1], w.ep[1])
  xmax = max(w.sp[1], w.ep[1])
  ymin = min(w.sp[2], w.ep[2])
  ymax = max(w.sp[2], w.ep[2])
  return xmin, ymin, xmax, ymax
end

function cellsize(a::Antidot)
  if a.where
    xmin = ymin = Inf
    xmax = ymax = -Inf
  else
    xmin, ymin = a.c .- a.r
    xmax, ymax = a.c .+ a.r
  end
  return xmin, ymin, xmax, ymax
end

"""
    cellsize(bt::Vector{Obstacle})
Return the delimiters `xmin, ymin, xmax, ymax` of the given billiard table.
"""
function cellsize(bt::Vector{Obstacle})

  xmin = ymin = Inf
  xmax = ymax = -Inf
  for obst in bt
    xs, ys, xm, ym = cellsize(obst)
    xmin = xmin > xs ? xs : xmin
    ymin = ymin > ys ? ys : ymin
    xmax = xmax < xm ? xm : xmax
    ymax = ymax < ym ? ym : ymax
  end
  return xmin, ymin, xmax, ymax
end


"""
```julia
randominside(bt::Vector{Obstacle}[, omega])
```
Return a particle with correct (allowed) initial conditions inside the given
billiard table defined by the vector `bt`. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `omega`.
Else, it is `Particle`.
"""
function randominside(bt::Vector{Obstacle})

  xmin, ymin, xmax, ymax = cellsize(bt)
  f = rand()
  while f == 0 || f==1/4 || f==1/2 || f == 3/4
    f = rand()
  end
  φ0 = f*2π

  xp = rand()*(xmax-xmin) + xmin
  yp = rand()*(ymax-ymin) + ymin
  p = Particle([xp, yp, φ0])

  dist = distance(p, bt)
  while dist <= 1e-12

    xp = rand()*(xmax-xmin) + xmin
    yp = rand()*(ymax-ymin) + ymin
    p.pos = [xp, yp]
    dist = distance(p, bt)
  end

  return p
end

function randominside(ω::Real, bt::Vector{Obstacle})

  if ω == 0
    return randominside(bt)
  end

  xmin, ymin, xmax, ymax = cellsize(bt)
  f = rand()
  while f == 0 || f==1/4 || f==1/2 || f == 3/4
    f = rand()
  end
  φ0 = f*2π

  xp = rand()*(xmax-xmin) + xmin
  yp = rand()*(ymax-ymin) + ymin
  p = MagneticParticle([xp, yp, φ0], ω)

  dist = distance(p, bt)
  while dist <= 1e-12

    xp = rand()*(xmax-xmin) + xmin
    yp = rand()*(ymax-ymin) + ymin
    p.pos = [xp, yp]
    dist = distance(p, bt)
  end

  return p
end
randominside(bt::Vector{Obstacle}, ω::Real) = randominside(ω, bt)
