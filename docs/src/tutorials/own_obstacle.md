# Creating your own `Obstacle` Type

In this tutorial we will go through the processes of creating a new obstacle type, a `Semicircle`. This type is already used in the [`billiard_bunimovich`](@ref) and [`billiard_mushroom`](@ref) functions.

!!! info "Everything uses `SVector{2}`"
    Fields of `Particle`s and `Obstacle`s contain all their information in 2-dimensional static vectors from module `StaticArrays`. This is important to keep in mind when extending new methods.

!!! info "Extends internal APIs"
    Notice that implementing your own obstacle requires you to extend methods that _do not_ belong to the public API.

!!! note "See also the `Ellipse` PR"
    Pull Request [#159](https://github.com/JuliaDynamics/DynamicalBilliards.jl/pull/159) implements the [`Ellipse`](@ref) obstacle, step by step, by following the tutorial of this page. All commits are commented and it can be a second helpful guide on how to implement an obstacle.

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
The following functions must obtain methods for `Semicircle` (or any other custom
`Obstacle`) in order for it to work with `DynamicalBilliards`:

1. [`DynamicalBilliards.normalvec`](@ref)
2. [`DynamicalBilliards.distance`](@ref) (with arguments `(position, obstacle)`)
3. [`DynamicalBilliards.collision`](@ref) with `Particle`

Assuming that upon collision a specular reflection happens, then you don't need
to define a method for [`DynamicalBilliards.specular!`](@ref). You can however define
custom methods for [`DynamicalBilliards.specular!`](@ref), which is what we have done e.g.
for [`RandomDisk`](@ref).

!!! note "Use `import`!"
    Notice that you have to properly `import` the methods to extend them. For example,
    do `import DynamicalBilliards: normalvec, distance, collision` time.

The first method is very simple, just do:
```julia
import DynamicalBilliards: normalvec, distance, collision
normalvec(d::Semicircle, pos) = normalize(d.c - pos)
```
Since the function is only used during `distance` and
[`DynamicalBilliards.resolvecollision!`](@ref) and since we will be writing explicit methods for the first,
we don't have to care about
what happens when the particle is far away from the boundary.

The `distance` method is a bit tricky. Since the type already subtypes `Circular`,
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

Finally, the method for [`collision`](@ref) is by far the most *trickiest*. But,
with pen, paper and a significant amount of patience, one can find a way:
```julia
function collision(p::Particle{T}, d::Semicircle{T})::T where {T}

    dc = p.pos - d.c
    B = dot(p.vel, dc)         #velocity towards circle center: B > 0
    C = dot(dc, dc) - d.r*d.r    #being outside of circle: C > 0
    Δ = B^2 - C

    Δ ≤ 0 && return nocollision(T)
    sqrtD = sqrt(Δ)

    nn = dot(dc, d.facedir)
    if nn ≥ 0 # I am NOT inside semicircle
        # Return most positive time
        t = -B + sqrtD
    else # I am inside semicircle:
        t = -B - sqrtD
        # these lines make sure that the code works for ANY starting position:
        if t ≤ 0 || distance(p, d) ≤ accuracy(T)
            t = -B + sqrtD
        end
    end
    # This check is necessary to not collide with the non-existing side
    newpos = p.pos + p.vel * t
    if dot(newpos - d.c, d.facedir) ≥ 0 # collision point on BAD HALF;
        return nocollision(T)
    end
    # If collision time is negative, return Inf:
    t ≤ 0.0 ? nocollision(T) : (t, p.pos + t*p.vel)
end
```

And that is all. The obstacle now works perfectly fine for straight propagation.

!!! note "Ray-Splitting support"
    Supporting ray-splitting for your custom obstacle is very easy. The first step is to give it a field called `pflag`, which is a `Bool`. The second step is to ensure that `collisiontime` works properly for particles coming from both directions of the obstacle! Both inside or outside! This is implemented for `Ellipse` in Pull Request [#159](https://github.com/JuliaDynamics/DynamicalBilliards.jl/pull/159).

## Optional Methods

1. [`DynamicalBilliards.cellsize`](@ref) : Enables [`randominside`](@ref) with this obstacle.
1. [`collision`](@ref) with [`MagneticParticle`](/basic/high_level/#particles) : enables magnetic propagation
2. [`plot`](@ref) with `obstacle` : enables plotting
3. [`DynamicalBilliards.specular!`](@ref) with `offset` : Allows [`lyapunovspectrum`](@ref) to be computed.
4. [`to_bcoords`](@ref) : Allows the [`boundarymap`](@ref) and [`boundarymap_portion`](@ref) to be computed.
5. [`from_bcoords`](@ref) : Allows [`phasespace_portion`](@ref) to be computed.

The [`DynamicalBilliards.cellsize`](@ref) method is kinda trivial:
```julia
import DynamicalBilliards: cellsize, plot, to_bcoords, from_bcoords

function cellsize(a::Semicircle{T}) where {T}
    xmin, ymin = a.c - a.r
    xmax, ymax = a.c + a.r
    return xmin, ymin, xmax, ymax
end
```


The [`collision`](@ref) method for [`MagneticParticle`](/basic/high_level/#particles) is also tricky, however it is almost identical with the method for the general `Circular` obstacle:
```julia
function collision(p::MagneticParticle{T}, o::Semicircle{T})::T where {T}
    ω = p.omega
    pc, rc = cyclotron(p)
    p1 = o.c
    r1 = o.r
    d = norm(p1-pc)
    if (d >= rc + r1) || (d <= abs(rc-r1))
        return nocollision(T)
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
    # Collision time, equiv. to arc-length until collision point:
    θ, I = nocollision(T)
    if cond1 || cond2
        for (Y, cond) in ((I1, cond1), (I2, cond2))
            if cond
                φ = realangle(p, o, Y)
                φ < θ && (θ = φ; I = Y)
            end
        end
    end
    # Collision time = arc-length until collision point
    return θ*rc, I
end

```

Then, we add swag by writing a method for [`plot`](@ref):

```julia
using PyPlot

function plot(d::Semicircle; kwargs...)
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

To enable computation of Lyapunov exponents in billiards with your obstacle,
you have to write another method for `specular!` that also handles the evolution
of perturbation vectors in tangent space. For this, the method has to
accept an argument of type `Vector{SVector{4, T}}`, which contains the four
perturbation vectors corresponding to the four Lyapunov exponents.

Finding a formula for the evolution of the perturbations requires some tricky
calculations. Fortunately for us, the results for general circular obstacles
were already determined by Dellago, Posch and Hoover [1] – we just have to
implement them.

```julia
function specular!(p::Particle{T}, o::Circular{T},
	offset::Vector{SVector{4, T}}) where {T<:AbstractFloat}

    n = normalvec(o, p.pos)
    ti = SV{T}(-p.vel[2],p.vel[1])

    cosa = -dot(n, p.vel)
    p.vel = p.vel + 2*cosa*n

    tf = SV{T}(-p.vel[2], p.vel[1])

    for k in 1:4
        δqprev = offset[k][δqind]
        δpprev = offset[k][δpind]
        # Formulas from Dellago, Posch and Hoover, PRE 53, 2, 1996: 1485-1501 (eq. 27)
        # with norm(p) = 1
        δq  = δqprev - 2*dot(δqprev,n)*n
        δp  = δpprev - 2*dot(δpprev,n)*n - curvature(o)*2/o.r*dot(δqprev,ti)/cosa*tf
        ###
        offset[k] = vcat(δq, δp)
    end
end

@inline curvature(::Semicircle) = -1
@inline curvature(::Disk) = +1

```

Note that calculating Lyapunov exponents for magnetic particles requires a
separate method, as the formulas are different for magnetic propagation.

Finally, we also add a methods for [`to_bcoords`](@ref) and [`from_bcoords`](@ref).
For them, see the relevant source file (use `@which`).


**References**

[1] : Ch. Dellago et al., Phys. Rev. E 53, 1485 (1996).
