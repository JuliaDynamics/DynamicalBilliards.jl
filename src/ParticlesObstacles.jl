using StaticArrays
import Base.show

export AbstractParticle, Particle, MagneticParticle, magnetic2standard,
standard2magnetic, cyclotron, Obstacle, Disk, Antidot, RandomDisk,
FiniteWall, PeriodicWall, RandomWall, SplitterWall,
normalvec, randominside, distance, cellsize, Wall, Circular

####################################################
## Particles
####################################################
"""
    AbstractParticle
Particle supertype.
"""
abstract type AbstractParticle{T<:AbstractFloat} end

"""
    Particle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table (numbers of type T).
#### Fields:
* `pos::SVector{2,T}` : Current position vector.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at.
  (Used only in periodic billiards)
#### Additional constructors:
```julia
Particle(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x, y, φ)
Particle() = Particle(rand(), rand(), rand()*2π)
```
"""
mutable struct Particle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    Particle{T}(pos, vel, cc) where {T<:AbstractFloat} = new(pos, normalize(vel), cc)
end

function Particle{T<:AbstractFloat}(ic::AbstractVector{T})
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
    return Particle(pos, vel, SVector{2,T}(0,0))
end
Particle(x, y, φ) = Particle(collect(promote(x,y,φ)))
Particle() = Particle(rand(), rand(), rand()*2π)
show(io::IO, p::Particle{T}) where {T} =
print(io, "Particle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)")

"""
    MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table with perpendicular magnetic field.
#### Fields:
* `pos::SVector{2,T}` : Current position vector.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at
  (Used only in periodic billiards).
* `omega::T` : Angular velocity (cyclic frequency) of rotational motion.
  Radius of rotation is `r=1/omega`.
#### Additional constructors:
```julia
MagneticParticle(ic::AbstractVector{T}, ω::Real) #where ic = [x0, y0, φ0]
MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)
```
"""
mutable struct MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    omega::T
    function MagneticParticle{T}(pos, vel, current_cell, ω) where {T<:AbstractFloat}
        if omega==0
            throw(ArgumentError("Angular velocity of magnetic particle cannot be 0."))
        end
        new(pos, normalize(vel), current_cell, ω)
    end
end

function MagneticParticle{T<:Real}(ic::AbstractVector{T}, ω::Real)
    φ0 = ic[3]
    S = T<:Integer ? Float64 : T
    pos = SVector{2,S}(ic[1:2]); vel = SVector{2,S}(cos(φ0), sin(φ0))
    return MagneticParticle(pos, vel, SVector{2,S}(0,0), convert(S,ω))
end
function MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
    a = collect(promote(x0, y0, φ0, ω))
    MagneticParticle(a[1:3], a[4])
end
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)

show(io::IO, p::MagneticParticle{T}) where {T} =
print(io, "MagneticParticle{$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)\nang. velocity: $(p.omega)")

"""
```julia
magnetic2standard(p::MagneticParticle, use_cell = true)
```
Create a standard `Particle` from a `MagneticParticle`.
"""
function magnetic2standard(p::MagneticParticle{T}, use_cell = true) where {T}
    pos = p.pos
    cell = use_cell ? current_cell : SVector{T,2}(0,0)
    return Particle(pos, p.vel, cell)
end

"""
```julia
standard2magnetic(p::Particle, omega, use_cell = true)
```
Create a `MagneticParticle` from a `Particle`.
"""
function standard2magnetic(p::Particle{T}, ω::Real, use_cell = true) where {T}
    pos = p.pos
    cell = use_cell ? current_cell : SVector{T,2}(0,0)
    return MagneticParticle(pos, p.vel, cell, convert(T, ω))
end


"""
```julia
cyclotron(p::MagneticParticle, use_cell = false)
```
Return center and radius of circular motion performed by the particle based on
`p.pos` (or `p.pos + p.current_cell`) and `p.vel`.
"""
function cyclotron(p::MagneticParticle{T}, use_cell = false) where {T}
    ω = p.omega
    pos = use_cell ? p.pos + p.current_cell : p.pos
    c::SVector{2, T} = pos - (1/ω)*[p.vel[2], -p.vel[1]]
    r = abs(1/ω)
    return c, r
end

####################################################
## Obstacles
####################################################
"""
    Obstacle{<:AbstractFloat}
Obstacle supertype.
"""
abstract type Obstacle{T<:AbstractFloat} end

"""
    Circular{T<:AbstractFloat} <: Obstacle{T}
Circular obstacle supertype.
"""
abstract type Circular{T<:AbstractFloat} <: Obstacle{T} end

"""
    Disk{T<:AbstractFloat}  <: Circular{T}
Disk-like obstacle with propagation allowed outside of the circle.
#### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `name::String` : Some name given for user convenience. Defaults to "Disk".
"""
struct Disk{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    name::String
end
function Disk(c::AbstractVector{T}, r::Real, name::String = "Disk") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Disk{S}(SVector{2,S}(c), convert(S, abs(r)), name)
end

"""
    RandomDisk{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that randomly (and uniformly) reflects colliding particles.
The propagation is allowed outside of the circle.
#### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `name::String` : Some name given for user convenience. Defaults to "Random disk".
"""
struct RandomDisk{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    name::String
end
function RandomDisk(
    c::AbstractVector{T}, r::Real, name::String = "RandomDisk") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return RandomDisk{S}(SVector{2,S}(c), convert(S, abs(r)), name)
end

"""
    Antidot{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that allows propagation both inside and outside of the disk.
Used in ray-splitting billiards.
#### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `pflag::Bool` : Flag that keeps track of where the particle is currently
  propagating (`pflag` = propagation-flag).
  `true` stands for *outside* the disk, `false` for *inside* the disk.
  Defaults to `true`.
* `name::String` : Name of the obstacle given for user convenience.
  Defaults to "Antidot".
"""
mutable struct Antidot{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    pflag::Bool
    name::String
end
function Antidot(c::AbstractVector{T}, r::Real,
    pflag = true, name::String = "Antidot") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Antidot{S}(SVector{2,S}(c), convert(S, abs(r)), pflag, name)
end

show(io::IO, w::Circular{T}) where {T} =
print(io, "$(w.name) ($T)\n", "center: $(w.c)\nradius: $(w.r)")


"""
    Wall{T<:AbstractFloat} <: Obstacle{T}
Wall obstacle supertype.
"""
abstract type Wall{T<:AbstractFloat} <: Obstacle{T} end

"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision.
#### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside of the
  billiard table). The size of the vector is irrelevant
  since it is internally normalized.
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Wall".
"""
struct FiniteWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    name::String
    function FiniteWall{T}(sp, ep, normal, name) where {T<:AbstractFloat}
        n = normalize(normal)
        d = dot(n, ep-sp)
        if abs(d) > 10eps(T)
            error("Normal vector is not actually normal to the wall")
        end
        return new(sp, ep, n, name)
    end
end
function FiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Wall")

    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return FiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

"""
    RandomWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing (uniformly) random reflection during collision.
#### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *is expected to come from* (pointing towards the inside of the
  billiard table).
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Random wall".
"""
struct RandomWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    name::String
    function RandomWall{T}(sp, ep, normal, name) where {T<:AbstractFloat}
        n = normalize(normal)
        d = dot(n, ep-sp)
        if abs(d) > 10eps(T)
            error("Normal vector is not actually normal to the wall")
        end
    return new(sp, ep, n, name)
  end
end
function RandomWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Random wall")

    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return RandomWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

"""
    PeriodicWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle that imposes periodic boundary conditions upon collision.
#### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from* (to the inside the billiard table).
  The size of the vector is **important**!
  This vector is added to a particle's `pos` during collision. Therefore the
  size of the normal vector must be correctly associated with the size of the
  periodic cell.
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Periodic wall".
"""
struct PeriodicWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    name::String
    function PeriodicWall{T}(sp, ep, normal, name) where {T<:AbstractFloat}
        d = dot(n, ep-sp)
        if abs(d) > 10eps(T)
            error("Normal vector is not actually normal to the wall")
        end
        return new(sp, ep, normal, name)
    end
end
function PeriodicWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Periodic wall")

    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return PeriodicWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end


"""
    SplitterWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing allowing for ray-splitting.
#### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision*.
  The size of the vector is irrelevant.
* `pflag::Bool` : Flag that keeps track of where the particle is currently
  propagating (`pflag` = propagation flag).
  `true` is associated with the `normal` vector the wall is
  instantiated with. Defaults to `true`.
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Splitter wall".
"""
mutable struct SplitterWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    pflag::Bool
    name::String
    function SplitterWall{T}(sp, ep, normal, pflag, name) where {T<:AbstractFloat}
        n = normalize(normal)
        d = dot(n, ep-sp)
        if abs(d) > 10eps(T)
            error("Normal vector is not actually normal to the wall")
        end
    return new(sp, ep, n, pflag, name)
  end
end
function SplitterWall(sp::AbstractVector, ep::AbstractVector,
    normal::AbstractVector, name::String = "Splitter wall")

    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return SplitterWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

#pretty print:
show(io::IO, w::Wall{T}) where {T} = print(io, "$(w.name) ($T)\n",
"start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")

"""
```julia
normalvec(obst::Obstacle, position)
```
Return the vector normal to the obstacle's boundary at the given position (which is
assumed to be very close to the obstacle's boundary).
"""
normalvec(wall::Wall, pos) = wall.normal
normalvec(w::PeriodicWall, pos) = normalize(w.normal)
normalvec(w::SplitterWall, pos) = w.pflag ? w.normal : -w.normal
normalvec(disk::Circular, pos) = normalize(pos - disk.c)
normalvec(a::Antidot, pos) = a.pflag ? normalize(pos - a.c) : -normalize(pos - a.c)

####################################################
## Distances
####################################################
"""
```julia
distance(p::AbstractParticle, o::Obstacle)
```
Return the **signed** distance between particle `p` and obstacle `o`, based on
`p.pos`. Positive distance corresponds to the particle being on the *allowed* region
of the `Obstacle`. E.g. for a `Disk`, the distance is positive when the particle is
outside of the disk, negative otherwise.

```julia
distance(p::AbstractParticle, bt::Vector{<:Obstacle})
```
Return minimum `distance(p, obst)` for all `obst` in `bt`.
If the `distance(p, bt)` is negative this means that the particle is outside
the billiard table.
"""
function distance(
    p::AbstractParticle{T}, bt::Vector{<:Obstacle{T}})::T where {T<:AbstractFloat}
    mindist::T = T(Inf)
    for i in eachindex(bt)
        dist::T = distance(p, bt[i])
        if dist <= mindist; mindist = dist; end
    end
    return mindist
end

function distance(p::AbstractParticle, w::Wall)
    v1 = p.pos - w.sp
    dot(v1, normalvec(w, p.pos))
end

# no new distance needed for SplitterWall because the `pflag` field
# has the necessary information to give the correct dinstance,
# since the distance is calculated throught the normalvec.

distance(p::AbstractParticle, d::Circular) = norm(p.pos - d.c) - d.r

function distance(p::AbstractParticle, a::Antidot)
    d = norm(p.pos - a.c) - a.r
    a.pflag ? d : -d
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

Used in `randominside()`, error checking and plotting.
"""
function cellsize(d::Circular{T}) where {T}
    xmin = ymin = T(Inf)
    xmax = ymax = T(-Inf)
    return xmin, ymin, xmax, ymax
end

function cellsize(w::Wall)
    xmin = min(w.sp[1], w.ep[1])
    xmax = max(w.sp[1], w.ep[1])
    ymin = min(w.sp[2], w.ep[2])
    ymax = max(w.sp[2], w.ep[2])
    return xmin, ymin, xmax, ymax
end

function cellsize(a::Antidot{T}) where {T}
    if a.pflag
        xmin = ymin = T(Inf)
        xmax = ymax = T(-Inf)
    else
        xmin, ymin = a.c .- a.r
        xmax, ymax = a.c .+ a.r
    end
    return xmin, ymin, xmax, ymax
end

function cellsize(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}

    xmin::T = ymin::T = T(Inf)
    xmax::T = ymax::T = T(-Inf)
    for i in eachindex(bt)
        xs::T, ys::T, xm::T, ym::T = cellsize(bt[i])
        xmin = xmin > xs ? xs : xmin
        ymin = ymin > ys ? ys : ymin
        xmax = xmax < xm ? xm : xmax
        ymax = ymax < ym ? ym : ymax
    end
    return xmin, ymin, xmax, ymax
end


"""
```julia
randominside(bt::Vector{<:Obstacle{T}}[, omega])
```
Return a particle with correct (allowed) initial conditions inside the given
billiard table defined by the vector `bt`. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `omega`.
Else, it is `Particle`.
"""
function randominside(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = f*2π

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = Particle([xp, yp, φ0])

    dist = distance(p, bt)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        p.pos = [xp, yp]
        dist = distance(p, bt)
    end

    return p
end

function randominside(ω::Real, bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    ω = convert(T, ω)
    if ω == 0
        return randominside(bt)
    end

    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = f*2π

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = MagneticParticle([xp, yp, φ0], ω)

    dist = distance(p, bt)
    while dist <= sqrt(eps(T))

        xp = rand()*(xmax-xmin) + xmin
        yp = rand()*(ymax-ymin) + ymin
        p.pos = [xp, yp]
        dist = distance(p, bt)
    end

    return p
end
randominside(bt::Vector{<:Obstacle{T}}, ω::Real) where {T} = randominside(ω, bt)
