using StaticArrays
import Base: show, eltype, getindex

export AbstractParticle, Particle, MagneticParticle,
cyclotron, Obstacle, Disk, Antidot, RandomDisk,
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
eltype(p::AbstractParticle{T}) where {T} = T

"""
    Obstacle{<:AbstractFloat}
Obstacle supertype.
"""
abstract type Obstacle{T<:AbstractFloat} end
eltype(o::Obstacle{T}) where {T} = T

eltype(bt::Vector{<:Obstacle{T}}) where {T} = T

"""
    Particle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table (mutable type).
### Fields:
* `pos::SVector{2,T}` : Current position vector.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at.
  (Used only in periodic billiards)
### Additional constructors:
```julia
Particle(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x, y, φ)
Particle() = Particle(rand(), rand(), rand()*2π)
Particle(bt::Vector{<:Obstacle}) = randominside(bt)
```
"""
mutable struct Particle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    function Particle(
        pos::SVector{2,T}, vel::SVector{2,T}, cc::SVector{2,T}) where{T<:AbstractFloat}
        new{T}(pos, normalize(vel), cc)
    end
end

function Particle(ic::AbstractVector{S}) where {S<:Real}
    T = S<:Integer ? Float64 : S
    φ0 = ic[3]
    pos = SVector{2,T}(ic[1:2]); vel = SVector{2,T}(cos(φ0), sin(φ0))
    return Particle(pos, vel, SVector{2,T}(0,0))
end
Particle(x::Real, y::Real, φ::Real) = Particle(collect(promote(x,y,φ)))
Particle() = Particle(rand(), rand(), rand()*2π)
Particle(bt::Vector{<:Obstacle}) = randominside(bt)
show(io::IO, p::Particle{T}) where {T} =
print(io, "Particle {$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)")

"""
    MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
Two-dimensional particle in a billiard table with perpendicular magnetic field
(mutable type).
### Fields:
* `pos::SVector{2,T}` : Current position vector.
* `vel::SVector{2,T}` : Current velocity vector (always of measure 1).
* `current_cell::SVector{2,T}` : Current "cell" the particle is located at
  (Used only in periodic billiards).
* `omega::T` : Angular velocity (cyclic frequency) of rotational motion.
  Radius of rotation is `r=1/omega`.
### Additional constructors:
```julia
MagneticParticle(ic::AbstractVector{T}, ω::Real) #where ic = [x0, y0, φ0]
MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)
MagneticParticle(bt::Vector{<:Obstacle}, ω) = randominside(bt, ω)
```
"""
mutable struct MagneticParticle{T<:AbstractFloat} <: AbstractParticle{T}
    pos::SVector{2,T}
    vel::SVector{2,T}
    current_cell::SVector{2,T}
    omega::T
    function MagneticParticle(pos::SVector{2,T}, vel::SVector{2,T},
        current_cell::SVector{2,T}, ω::T) where {T<:AbstractFloat}
        if ω==0
            throw(ArgumentError("Angular velocity of magnetic particle cannot be 0."))
        end
        new{T}(pos, normalize(vel), current_cell, ω)
    end
end

function MagneticParticle(ic::AbstractVector{T}, ω::Real) where {T<:Real}
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
MagneticParticle(bt::Vector{<:Obstacle}, ω) = randominside(bt, ω)

show(io::IO, p::MagneticParticle{T}) where {T} =
print(io, "Magnetic particle {$T}\n",
"position: $(p.pos+p.current_cell)\nvelocity: $(p.vel)\nang. velocity: $(p.omega)")


"""
    cyclotron(p::MagneticParticle, use_cell = false)
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
    Circular{T<:AbstractFloat} <: Obstacle{T}
Circular obstacle supertype.
"""
abstract type Circular{T<:AbstractFloat} <: Obstacle{T} end

"""
    Disk{T<:AbstractFloat}  <: Circular{T}
Disk-like obstacle with propagation allowed outside of the circle (immutable type).
### Fields:
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
Disk{T}(args...) where {T} = Disk(args...)

"""
    RandomDisk{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that randomly (and uniformly) reflects colliding particles.
The propagation is allowed outside of the circle.
### Fields:
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
RandomDisk{T}(args...) where {T} = RandomDisk(args...)

"""
    Antidot{T<:AbstractFloat} <: Circular{T}
Disk-like obstacle that allows propagation both inside and outside of the disk
(immutable type). Used in ray-splitting billiards.
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
    pflag::Bool = true, name::String = "Antidot") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Antidot{S}(SVector{2,S}(c), convert(S, abs(r)), pflag, name)
end
Antidot(c, r, name::String = "Antidot") = Antidot(c, r, true, name)
Antidot{T}(args...) where {T} = Antidot(args...)

show(io::IO, w::Circular{T}) where {T} =
print(io, "$(w.name) {$T}\n", "center: $(w.c)\nradius: $(w.r)")


"""
    Wall{T<:AbstractFloat} <: Obstacle{T}
Wall obstacle supertype.
"""
abstract type Wall{T<:AbstractFloat} <: Obstacle{T} end

"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision (immutable type).
### Fields:
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
end
function FiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Wall")
    T = eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return FiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

"""
    RandomWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing (uniformly) random reflection during collision (immutable type).
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
end
function RandomWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Random wall")
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    return RandomWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

"""
    PeriodicWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle that imposes periodic boundary conditions upon collision (immutable type).
### Fields:
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
end
function PeriodicWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Periodic wall")
    T = eltype(sp)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return PeriodicWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end


"""
    SplitterWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing allowing for ray-splitting (immutable type).
### Fields:
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
end
function SplitterWall(sp::AbstractVector, ep::AbstractVector,
    normal::AbstractVector, pflag::Bool = true, name::String = "Splitter wall")
    T = eltype(sp)
    n = normalize(normal)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return SplitterWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), pflag, name)
end
SplitterWall(sp, ep, n, name::String = "Splitter wall") = SplitterWall(sp, ep, n, true, name)
#pretty print:
show(io::IO, w::Wall{T}) where {T} = print(io, "$(w.name) {$T}\n",
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
## Billiard Table
####################################################
function isperiodic(bt)::Bool
    for obst in bt
        if typeof(obst) <: PeriodicWall
            return true
        end
    end
    return false
end

immutable BilliardTable{T, BT<:Tuple}
    bt::BT
end

function BilliardTable(bt)
    T = eltype(bt[1])
    if typeof(bt) <: Tuple
        return BilliardTable{T, typeof(bt)}(bt)
    elseif typeof(bt) <: AbstractVector
        tup = (bt...)
        return BilliardTable{T, typeof(tup)}(tup)
    else
        throw(ArgumentError("Argument given to BilliardTable is not a container"))
    end
end

getobstacle(bt::BilliardTable{T,S}, ::Val{N}) where {T,S,N} = bt.bt[N]


ObstInd(::Val{N}) where {N} = SVector{1, Val{N}}(Val(N))
#This doesn't work:
getobstacle(bt::BilliardTable{T,S}, N::Int) where {T,S} =
getobstacle(bt, ObstInd(Val(N)))

getindex(bt::BilliardTable, i) = bt.bt[i]

####################################################
## Distances
####################################################
"""
    distance(p::AbstractParticle, o::Obstacle)
Return the **signed** distance between particle `p` and obstacle `o`, based on
`p.pos`. Positive distance corresponds to the particle being on the *allowed* region
of the `Obstacle`. E.g. for a `Disk`, the distance is positive when the particle is
outside of the disk, negative otherwise.

    distance(p::AbstractParticle, bt::Vector{<:Obstacle})
Return minimum `distance(p, obst)` for all `obst` in `bt`.
If the `distance(p, bt)` is negative this means that the particle is outside
the billiard table.

All `distance` functions can also be given a position (Vector) instead of a particle.
"""
function distance(pos::AbstractVector{T}, v::Vector{<:Obstacle{T}})::T where {T}
    d = T(Inf)
    for obst in v
        di = distance(pos, obst)
        di < d && (d = di)
    end
    return d
end

(distance(p::AbstractParticle{T}, v::Vector{<:Obstacle{T}})::T) where {T} =
distance(p.pos, v)

(distance(pos::AbstractVector{T}, bt::Tuple)::T) where {T<:AbstractFloat} =
min(distance(pos, obst) for obst in bt)

distance(pos::AbstractVector, bt::BilliardTable) = distance(pos, bt.bt)
distance(p::AbstractParticle, bt::BilliardTable) = distance(p.pos, bt.bt)

function distance(pos::AbstractVector{T}, w::Wall{T})::T where {T}
    v1 = pos - w.sp
    dot(v1, normalvec(w, pos))
end

# no new distance needed for SplitterWall because the `pflag` field
# has the necessary information to give the correct dinstance,
# since the distance is calculated through the normalvec.

distance(pos::AbstractVector{T}, d::Circular{T}) where {T} = norm(pos - d.c) - d.r

function distance(
    pos::AbstractVector{T}, a::Antidot{T}, pflag::Bool)::T where {T}
    d = norm(pos - a.c) - a.r
    a.pflag ? d : -d
end
distance(pos::AbstractVector{T}, a::Antidot{T}) where {T} = distance(pos, a, a.pflag)

(distance(p::AbstractParticle{T}, obst::Obstacle{T})::T) where {T} =
distance(p.pos, obst)


####################################################
## Initial Conditions
####################################################

"""
    cellsize(bt)
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
    randominside(bt::Vector{<:Obstacle{T}}[, ω])
Return a particle with correct (allowed) initial conditions inside the given
billiard table defined by the vector `bt`. If supplied with a second argument the
type of the returned particle is `MagneticParticle`, with angular velocity `ω` (unless
`ω` is 0). Else, it is `Particle`.
"""
function randominside(bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}
    xmin::T, ymin::T, xmax::T, ymax::T = cellsize(bt)
    f = T(rand())
    while f == 0 || f==1/4 || f==1/2 || f == 3/4
        f = T(rand())
    end
    φ0 = T(f*2π)

    xp = T(rand())*(xmax-xmin) + xmin
    yp = T(rand())*(ymax-ymin) + ymin
    p = Particle([xp, yp, φ0])

    dist = distance(p, bt)
    while dist <= sqrt(eps(T))

        xp = T(rand())*(xmax-xmin) + xmin
        yp = T(rand())*(ymax-ymin) + ymin
        p.pos = SVector{2,T}(xp, yp)
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
    φ0 = T(f*2π)

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
