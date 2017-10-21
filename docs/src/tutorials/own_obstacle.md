# Creating your own `Obstacle` Type

In this tutorial we will go through the processes of creating a new obstacle type, a
`Semicircle`. This type is already used in the [`billiard_bunimovich`](@ref) and
[`billiard_mushroom`](@ref) functions.

## Type Definition
The first thing you have to do is make your new type a sub-type of `Obstacle{T}`
(or any other abstract sub-type of it). We will do:
```julia
struct Semicircle{T<:AbstractFloat} <: Circular{T} # <: Obstacle{T}
    c::SVector{2,T}
    r::T
    facedir::SVector{2,T}
    name::String
end
```
`c` is the center and `r` is the radius of the full circle. `facedir` is the direction
which the semicircle is facing, which is also the direction of its "open" face.
`name` is a final (mandatory) field that is there to easily identify which obstacle
is which.

Notice that the `struct` must be parameterized by `T<:AbstractFloat` (see
the [numerical precision](/physics/#numerical-precision) page for more).

For convenience, I will also define:
```julia
function Semicircle(
    c::AbstractVector{T}, r::Real, facedir, name = "Semicircle") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Semicircle{S}(SVector{2,S}(c), convert(S, abs(r)), name)
end
```
so that constructing a `Semicircle` is possible from arbitrary vectors.

## Necessary Methods
The following functions must obtain methods for `Semicircle` (or any other custom
`Obstacle`) in order for it to
work with `DynamicalBilliards.jl`:

1. [`normalvec`](@ref)
2. [`distance`](@ref) (with arguments `(position, obstacle)`)
3. [`collisiontime`](@ref) with `Particle`
2. Trivial `collisiontime` method with `MagneticParticle` that returns `error()`
   (see optional methods)

Assuming that upon collision a specular reflection happens, then you don't need
to define a method for [`resolvecollision!`](@ref).

The first method is very simple, just do:
```julia
normalvec(d::Semicircle, pos) = normalize(d.c - pos)
```
Since the function is only used during collision time estimation, `distance` and
`resolvecollision!`, and since we will be writing explicit methods for the first two,
we don't have to care about
what happens when the particle is far away from the boundary.

The `distance` method is a bit tricky. Since the type already subtypes `Circular`,
the following definition from `DynamicalBilliards` applies:
```julia
distance(pos::AbstractVector{T}, d::Circular{T}) where {T} = norm(pos - d.c) - d.r
```
However, the method must be
expanded. That is because when the particle is on the "open" half of the
disk, the distance is not correct. We write:
```julia
SV = SVector{2} #convenience
function distance(pos::AbstractVector{T}, s::Semicircle{T}) where {T}
    # Check on which half of circle is the particle
    v1 = pos .- s.c
    nn = dot(v1, s.facedir)
    if nn ≤ 0 # I am "inside semicircle
        return s.r - norm(pos - s.c)
    else # I am on the "other side"
        end1 = SV(s.c[1] + s.r*s.facedir[2], s.c[2] - s.r*s.facedir[1])
        end2 = SV(s.c[1] - s.r*s.facedir[2], s.c[2] + s.r*s.facedir[1])
        return min(norm(pos - end1), norm(pos - end2))
    end
end
```
Notice that this definition always returns positive distance when the particle is on
the "other side".

Finally, the method for [`collisiontime`](@ref) is by far the most *trickiest*. But,
with pen, paper and some patience, one can find a way:
```julia
function collisiontime(p::Particle{T}, d::Semicircle{T})::T where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         
    C = dot(dc, dc) - d.r^2    
    Δ = B^2 - C

    Δ <= 0 && return Inf
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        # these lines make sure that the code works for ANY starting position:
        t = -B - sqrtD
        if t ≤ 0
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel .* t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return Inf
    end
    # If collision time is negative, return Inf:
    t <= 0.0 ? Inf : t
end
```

And that is all. The obstacle now works perfectly fine for straight propagation
and properly initializes particles with `randominside`!



## Optional Methods

1. [`collisiontime`](@ref) with `MagneticParticle` : enables magnetic propagation
2. `plot_obstacle` : enables plotting (used in `plot_billiard`)


for 2:
```julia
Arc = PyPlot.matplotlib[:patches][:Arc]
function plot_obstacle(d::Semicircle; kwargs...)
  theta1 = atan2(d.facedir[2], d.facedir[1])*180/π + 90
  theta2 = theta1 + 180
  s1 = Arc(d.c, 2d.r, 2d.r, theta1 = theta1, theta2 = theta2,
  edgecolor = (0,0.6,0), linewidth = 2.0, kwargs...)
  PyPlot.gca()[:add_artist](s1)
  PyPlot.show()
end
```
