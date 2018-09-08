# Creating your own `Obstacle` Type

In this tutorial we will go through the processes of creating a new obstacle type, a
`Semicircle`. This type is already used in the [`billiard_bunimovich`](@ref) and
[`billiard_mushroom`](@ref) functions.

!!! info "Everything uses `SVector{2}`"
    Fields of `Particle`s and `Obstacle`s contain all their information in 2-dimensional static vectors from module `StaticArrays`. This is important to keep in mind when extending new methods.

## Type Definition
The first thing you have to do is make your new type a sub-type of `Obstacle{T}`
(or any other abstract sub-type of it). We will do:
```julia
struct Semicircle{T<:AbstractFloat} <: Circular{T} # <: Obstacle{T}
    c::SVector{2,T} # this MUST be a static vector
    r::T
    facedir::SVector{2,T} # this MUST be a static vector
    name::String # this is an OPTIONAL field
end
```
`c` is the center and `r` is the radius of the full circle. `facedir` is the direction
which the semicircle is facing, which is also the direction of its "open" face.

`name` is an optional field, that allows one to easily identify which obstacle
is which. It is also used when printing a [`Billiard`](@ref). If not used,
then `string(typeof(obstacle))` is used instead.

Notice that the `struct` **must be** parameterized by `T<:AbstractFloat` (see
the [numerical precision](/physics/#numerical-precision) page for more).

For convenience, we will also define:
```julia
function Semicircle(
    c::AbstractVector{T}, r::Real, facedir, name = "Semicircle") where {T<:Real}
    S = T <: Integer ? Float64 : T
    return Semicircle{S}(SVector{2,S}(c), convert(S, abs(r)), name)
end
```
so that constructing a `Semicircle` is possible from arbitrary vectors.

## Necessary Methods
The following functions must obtain methods for your custom
`Obstacle` in order for it to work with `DynamicalBilliards`:

1. [`normalvec`](@ref)
2. [`distance`](@ref) (with arguments `(position, obstacle)`)
3. [`collisiontime`](@ref) with `Particle`

Assuming that upon collision a specular reflection happens, then you don't need
to define a method for [`resolvecollision!`](@ref). You can however define
custom methods for [`resolvecollision!`](@ref), which is what we have done e.g.
for [`RandomDisk`](@ref).

!!! note "Use `import`!"
    Notice that you have to properly `import` the methods to extend them. For example,
    do `import DynamicalBilliards: normalvec, distance, collision` time.

The first method is very simple, just do:
```julia
import DynamicalBilliards: normalvec, distance, collisiontime
normalvec(d::Semicircle, pos) = normalize(d.c - pos)
```
Since the function is only used during [`distance`](@ref) and
[`resolvecollision!`](@ref) and since we will be writing explicit methods for the first,
we don't have to care about
what happens when the particle is far away from the boundary.

The [`distance`](@ref) method is a bit tricky. Since the type already subtypes `Circular`,
the following definition from `DynamicalBilliards` applies:
```julia
distance(pos::AbstractVector, d::Circular) = norm(pos - d.c) - d.r
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
    if nn ≤ 0 # I am "inside semicircle"
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
with pen, paper and a significant amount of patience, one can find a way:
```julia
function collisiontime(p::Particle{T}, d::Semicircle{T})::T where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r^2    #being outside of circle: C > 0
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
        if t ≤ 0 || distance(p, d) ≤ distancecheck(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel .* t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return Inf
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? Inf : t
end
```

And that is all. The obstacle now works perfectly fine for straight propagation.



## Optional Methods

1. [`cellsize`](@ref) : Enables [`randominside`](@ref) with this obstacle.
1. [`collisiontime`](@ref) with [`MagneticParticle`](/basic/high_level/#particles) : enables magnetic propagation
2. [`plot_obstacle`](@ref) : enables plotting (used in [`plot_billiard`](@ref)) (but requires [`cellsize`](@ref) to be already implemented, because [`plot_billiard`](@ref) also does automatic axis limits configuration)
1. [`translate`](@ref) : Enables plotting the obstacle with periodic billiards
3. [`to_bcoords`](@ref) : Allows the [`boundarymap`](@ref) and [`boundarymap_portion`](@ref) to be computed.
4. [`from_bcoords`](@ref) : Allows [`phasespace_portion`](@ref) to be computed.

The [`cellsize`](@ref) method is kinda trivial:
```julia
import DynamicalBilliards: cellsize, plot_obstacle, to_bcoords, from_bcoords

function cellsize(a::Semicircle{T}) where {T}
    xmin, ymin = a.c - a.r
    xmax, ymax = a.c + a.r
    return xmin, ymin, xmax, ymax
end
```


The [`collisiontime`](@ref) method for [`MagneticParticle`](/basic/high_level/#particles) is also
tricky, however it is almost identical with the method for the general [`Circular`](@ref) obstacle:
```julia
function collisiontime(p::MagneticParticle{T}, o::Semicircle{T})::T where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d >= rc + r1) || (d <= abs(rc-r1))
        return Inf
    end
    # Solve quadratic:
    a = (rc^2 - r1^2 + d^2)/2d
    h = sqrt(rc^2 - a^2)
    # Collision points (always 2):
    I1 = SVector{2, T}(
        pc[1] + a*(p1[1] - pc[1])/d + h*(p1[2] - pc[2])/d,
        pc[2] + a*(p1[2] - pc[2])/d - h*(p1[1] - pc[1])/d
    )
    I2 = SVector{2, T}(
        pc[1] + a*(p1[1] - pc[1])/d - h*(p1[2] - pc[2])/d,
        pc[2] + a*(p1[2] - pc[2])/d + h*(p1[1] - pc[1])/d
    )
    # Only consider intersections on the "correct" side of Semicircle:
    cond1 = dot(I1-o.c, o.facedir) < 0
    cond2 = dot(I2-o.c, o.facedir) < 0
    if cond1 || cond2
        # Calculate real angle until intersection:
        θ1 = cond1 ? realangle(p, o, I1) : T(Inf)
        θ2 = cond2 ? realangle(p, o, I2) : T(Inf)
        # Collision time, equiv. to arc-length until collision point:
        return min(θ1, θ2)*rc
    else
        return Inf
    end
end

```

Then, we add swag by writing a method for [`plot_obstacle`](@ref):

```julia
using PyPlot

function plot_obstacle(d::Semicircle; kwargs...)
    theta1 = atan(d.facedir[2], d.facedir[1])*180/π + 90
    theta2 = theta1 + 180
    edgecolor = DynamicalBilliards.obcolor(d)
    s1 = PyPlot.matplotlib[:patches][:Arc](
        d.c, 2d.r, 2d.r, theta1 = theta1, theta2 = theta2, edgecolor = edgecolor,
        lw = 2.0, kwargs...)
    PyPlot.gca()[:add_artist](s1)
    PyPlot.show()
end
```

Finally, we also add a methods for [`to_bcoords`](@ref) and [`from_bcoords`](@ref).
For them, see the [relevant source file](https://github.com/JuliaDynamics/DynamicalBilliards.jl/blob/master/src/boundary/boundarymap.jl).
