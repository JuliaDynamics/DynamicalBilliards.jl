export Obstacle, Disk, Antidot, RandomDisk, Wall, Circular,
InfiniteWall, PeriodicWall, RandomWall, SplitterWall, FiniteWall,
normalvec, distance, cellsize, Semicircle,
arclength, totallength, real_pos, real_coordinates
export translate
#######################################################################################
## Circles
#######################################################################################
"""
    Obstacle{<:AbstractFloat}
Obstacle supertype.
"""
abstract type Obstacle{T<:AbstractFloat} end
eltype(o::Obstacle{T}) where {T} = T


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
(mutable type). Used in ray-splitting billiards.
### Fields:
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

"""
    Semicircle{T<:AbstractFloat} <: Circular{T}
Obstacle that represents half a circle.
### Fields:
* `c::SVector{2,T}` : Center.
* `r::T` : Radius.
* `facedir::SVector{2,T}` : Direction where the open face of the Semicircle is facing.
* `name::String` : Name of the obstacle given for user convenience.
  Defaults to "Semicircle".
"""
struct Semicircle{T<:AbstractFloat} <: Circular{T}
    c::SVector{2,T}
    r::T
    facedir::SVector{2,T}
    name::String
end
function Semicircle(
    c::AbstractVector{T}, r::Real, facedir, name = "Semicircle") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Semicircle{S}(
    SV(c), convert(S, abs(r)), SV(normalize(facedir)), name)
end

show(io::IO, w::Circular{T}) where {T} =
print(io, "$(w.name) {$T}\n", "center: $(w.c)\nradius: $(w.r)")

show(io::IO, w::Semicircle{T}) where {T} =
print(io, "$(w.name) {$T}\n", "center: $(w.c)\nradius: $(w.r)\nfacedir: $(w.facedir)")



#######################################################################################
## Walls
#######################################################################################
"""
    Wall{T<:AbstractFloat} <: Obstacle{T}
Wall obstacle supertype.
"""
abstract type Wall{T<:AbstractFloat} <: Obstacle{T} end

"""
    InfiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision (immutable type).
Faster than [`FiniteWall`](@ref), meant to be used for convex billiards.
### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside of the
  billiard). The size of the vector is irrelevant
  since it is internally normalized.
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Wall".
"""
struct InfiniteWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    name::String
end
function InfiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, name::String = "Wall")
    T = eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    return InfiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), name)
end

"""
    FiniteWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing specular reflection during collision (immutable type).
Slower than [`InfiniteWall`](@ref), meant to be used for non-convex billiards.

Giving a `true` value to the field `isdoor` designates this obstacle to be a `Door`.
This is used in [`escapetime`](@ref) function. A `Door` is a obstacle of the
billiard that the particle can escape from, thus enabling calculations
of escape times.

### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *will come from before a collision* (pointing towards the inside of the
  billiard). The size of the vector is irrelevant
  since it is internally normalized.
* `isdoor::Bool` : Flag of whether this `FiniteWall` instance is a "Door".
* `name::String` : Name of the obstacle, given for user convenience.
  Defaults to "Finite Wall".
"""
struct FiniteWall{T<:AbstractFloat} <: Wall{T}
    sp::SVector{2,T}
    ep::SVector{2,T}
    normal::SVector{2,T}
    width::T
    center::SVector{2,T}
    isdoor::Bool
    name::String
end
function FiniteWall(sp::AbstractVector, ep::AbstractVector,
    n::AbstractVector, isdoor::Bool = false, name::String = "Finite Wall")
    T = eltype(sp)
    n = normalize(n)
    d = dot(n, ep-sp)
    if abs(d) > 10eps(T)
        error("Normal vector is not actually normal to the wall: dot = $d")
    end
    T = eltype(sp) <: Integer ? Float64 : eltype(sp)
    w = norm(ep - sp)
    center = @. (ep+sp)/2
    return FiniteWall{T}(SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n),
    w, SVector{2,T}(center), isdoor, name)
end

isdoor(w) = w.isdoor

"""
    RandomWall{T<:AbstractFloat} <: Wall{T}
Wall obstacle imposing (uniformly) random reflection during collision (immutable type).
### Fields:
* `sp::SVector{2,T}` : Starting point of the Wall.
* `ep::SVector{2,T}` : Ending point of the Wall.
* `normal::SVector{2,T}` : Normal vector to the wall, pointing to where the
  particle *is expected to come from* (pointing towards the inside of the
  billiard).
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
  particle *will come from* (to the inside the billiard).
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
Wall obstacle imposing allowing for ray-splitting (mutable type).
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
    return SplitterWall{T}(
    SVector{2,T}(sp), SVector{2,T}(ep), SVector{2,T}(n), pflag, name)
end
SplitterWall(sp, ep, n, name::String = "Splitter wall") =
SplitterWall(sp, ep, n, true, name)
#pretty print:
show(io::IO, w::Wall{T}) where {T} = print(io, "$(w.name) {$T}\n",
"start point: $(w.sp)\nend point: $(w.ep)\nnormal vector: $(w.normal)")


#######################################################################################
## Normal vectors
#######################################################################################
"""
    normalvec(obst::Obstacle, position)
Return the vector normal to the obstacle's boundary at the given position (which is
assumed to be very close to the obstacle's boundary).
"""
@inline normalvec(wall::Wall, pos) = wall.normal
@inline normalvec(w::PeriodicWall, pos) = normalize(w.normal)
@inline normalvec(w::SplitterWall, pos) = w.pflag ? w.normal : -w.normal
@inline normalvec(disk::Circular, pos) = normalize(pos - disk.c)
@inline normalvec(a::Antidot, pos) =
    a.pflag ? normalize(pos - a.c) : -normalize(pos - a.c)
@inline normalvec(d::Semicircle, pos) = normalize(d.c - pos)

#######################################################################################
## Arclengths
#######################################################################################

"""
    arclength(p::AbstractParticle, o::Obstacle)
Returns the boundary coordinate of the particle on the obstacle,
assuming that the particle position is on the obstacle.

The boundary coordinate is measured as:
* the distance from start point to end point in `Wall`s
* the arc length measured counterclockwise from the open face in `Semicircle`s
* the arc length measured counterclockwise from the rightmost point in `Circular`s
"""
arclength(p::AbstractParticle, o) = arclength(p.pos, o)
arclength(pos::SV, o::Wall) = norm(pos - o.sp)

function arclength(pos::SV{T}, o::Circular{T}) where {T<:AbstractFloat}
    #projecting pos to x Axis
    d = (pos - o.c)/o.r
    r = acos(clamp(d[1], -1, 1))*o.r
    return r
end

function arclength(pos::SV{T}, o::Semicircle{T}) where {T<:AbstractFloat}
    #project pos on open face
    chrd = SV(-o.facedir[2],o.facedir[1]) #tangent to open face
    d = (pos - o.c)/o.r
    x = dot(d, chrd)
    r =  acos(clamp(x, -1, 1))*o.r
    return r
end

#######################################################################################
## real_coordinates
#######################################################################################

function real_coordinates(ξ, sφ, o::Obstacle{T}) where T
    pos = real_pos(ξ, o)
    cφ = cos(asin(sφ))
    n = normalvec(o, pos)
    vel = SV{T}(-n[1]*cφ + n[2]*sφ, -n[1]*sφ - n[2]*cφ)
    return pos,vel
end


real_pos(ξ, o::Wall) = o.sp + ξ*normalize(o.ep - o.sp)

#TODO: Use `sincos` for julia 0.7
real_pos(ξ, o::Circular{T}) where T = o.c .+ o.r * SV{T}(cos(ξ/o.r), sin(ξ/o.r))

function real_pos(ξ, o::Semicircle{T}) where T
    θshift = acos(-o.facedir[2])
    return o.c .+ SV{T}(cos(ξ/o.r + θshift), sin(ξ/o.r + θshift))
end
"""
    totallength(o::Obstacle)
Return the total length of `o`.
"""
@inline totallength(o::Wall) = norm(o.ep - o.sp)
@inline totallength(o::Semicircle) = π*o.r
@inline totallength(o::Circular) = 2π*o.r



#######################################################################################
## Distances
#######################################################################################
"""
    project_to_line(point, c, n)
Project given `point` to line that contains point `c` and has **normal vector** `n`.
"""
@inline function project_to_line(point, c, n)
    posdot = dot(c - point, n)
    intersection = point + posdot*n
end

"""
    distance(p::AbstractParticle, o::Obstacle)
Return the **signed** distance between particle `p` and obstacle `o`, based on
`p.pos`. Positive distance corresponds to the particle being on the *allowed* region
of the `Obstacle`. E.g. for a `Disk`, the distance is positive when the particle is
outside of the disk, negative otherwise.

    distance(p::AbstractParticle, bd::Billiard)
Return minimum `distance(p, obst)` for all `obst` in `bd`.
If the `distance(p, bd)` is negative this means that the particle is outside
the billiard.

All `distance` functions can also be given a position (vector) instead of a particle.
"""
(distance(p::AbstractParticle{T}, obst::Obstacle{T})::T) where {T} =
distance(p.pos, obst)

@inline function distance(pos::AbstractVector{T}, w::Wall{T})::T where {T}
    v1 = pos - w.sp
    dot(v1, normalvec(w, pos))
end

# no new distance needed for SplitterWall because the `pflag` field
# has the necessary information to give the correct dinstance,
# since the distance is calculated through the normalvec.

@inline distance(pos::AbstractVector{T}, d::Circular{T}) where {T} =
    norm(pos - d.c) - d.r

@inline function distance(
    pos::AbstractVector{T}, a::Antidot{T}, pflag::Bool)::T where {T}
    d = norm(pos - a.c) - a.r
    a.pflag ? d : -d
end
@inline distance(pos::AbstractVector{T}, a::Antidot{T}) where {T} =
    distance(pos, a, a.pflag)

function distance(pos::AbstractVector{T}, s::Semicircle{T}) where {T}
    # Check on which half of circle is the particle
    v1 = pos - s.c
    nn = dot(v1, s.facedir)
    if nn ≤ 0 # I am "inside" semicircle
        return s.r - norm(pos - s.c)
    else # I am on the "other side"
        end1 = SV(s.c[1] + s.r*s.facedir[2], s.c[2] - s.r*s.facedir[1])
        end2 = SV(s.c[1] - s.r*s.facedir[2], s.c[2] + s.r*s.facedir[1])
        return min(norm(pos - end1), norm(pos - end2))
    end
end


# The entire functionality of `distance_init` is necessary only for
# FiniteWall !!!
distance_init(p::AbstractParticle, a::Obstacle) = distance_init(p.pos, a)
distance_init(pos::SVector, a::Obstacle) = distance(pos, a)

function distance_init(pos::SVector{2,T}, w::FiniteWall{T})::T where {T}

    n = normalvec(w, pos)
    posdot = dot(w.sp - pos, n)
    if posdot ≥ 0 # I am behind wall
        intersection = project_to_line(pos, w.center, n)
        dfc = norm(intersection - w.center)
        if dfc > w.width/2
            return +1.0 # but not directly behind
        else
            return -1.0
        end
    end
    v1 = pos - w.sp
    dot(v1, n)
end

####################################################
## Initial Conditions
####################################################
"""
    cellsize(bd)
Return the delimiters `xmin, ymin, xmax, ymax` of the given obstacle/billiard.

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

function cellsize(a::Semicircle{T}) where {T}
    xmin, ymin = a.c .- a.r
    xmax, ymax = a.c .+ a.r
    return xmin, ymin, xmax, ymax
end


####################################################
## Translations
####################################################
"""
    translate(obst::Obstacle, vector)
Create a copy of the given obstacle with its position
translated by `vector`.
"""
function translate end

for T in subtypes(Circular)
  @eval translate(d::$T, vec) = ($T)(d.c .+ vec, d.r)
end

for T in subtypes(Wall)
  @eval translate(w::$T, vec) = ($T)(w.sp + vec, w.ep + vec, w.normal)
end
