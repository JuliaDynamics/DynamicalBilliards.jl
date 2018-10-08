# Creating your own Billiard
```@docs
Billiard
```
---

A `Billiard` is a wrapper of a `Tuple` of `Obstacle`s.
The abstract Type `Obstacle{T}` is the supertype of all objects that a particle may collide with, with global billiard precision of type `T`.

There are many premade functions that construct well-known billiards, like the periodic Sinai billiard.
You can find all of them at the [Standard Billiards Library](/basic/high_level/#standard-billiards-library).

To create a custom billiard from scratch, it is often convenient to start with an empty `Vector{Obstacle{T}}`:
```@example tut1
using DynamicalBilliards
bd = Obstacle{Float64}[]  # T<: AbstractFloat
```
and then you create your obstacles one by one and add them to it. All obstacles that are already defined in the package
can be found at the [Obstacles page](/basic/library/#obstacles) of the library. The function [`billiard_polygon`](@ref) creates a polygonal billiard table.
However, for the example of this page, we will create a hexagonal billiard with a disk in the middle step-by-step.

The first step is to define the six walls of the billiard table.
An [`InfiniteWall`](@ref) object needs to be supplemented with a start point, an end point, a normal vector and, optionally, a name.

The vertex points of a regular hexagon of radius `r` are given by the formula:
```math
(x,y) = \left( r\cos\left(\frac{2\pi i}{6}\right), r\cos\left(\frac{2\pi i}{6}\right) \right)\,, \quad \text{for i $\in$ \{1,...,6\}}
```
To create each wall object, we will implement the following loop:
```@example tut1
hexagon_vertex = (r) -> [ [r*cos(2π*i/6), r*sin(2π*i/6)] for i in 1:6]
hexver = hexagon_vertex(2.0)

for i in eachindex(hexver)
  starting = hexver[i]
  ending = hexver[mod1(i+1, length(hexver))]
  w = ending - starting
  normal = [-w[2], w[1]]
  wall = InfiniteWall(starting, ending, normal, "wall $i")
  push!(bd, wall)
end

summary(bd)
```

The `normal` vector of a `Wall` obstacle is necessary to be supplemented by the user because it must point towards where the particle is expected to come from. If `w` is the vector (wall) pointing from start- to end-point then the vector `[-w[2], w[1]]` is pointing to the left of `w` and the vector `[w[2], -[w1]]` is pointing to the right. Both are normal to `w`, but you have to know which one to pick. In this case this is very easy, since the normal has to simply point towards the origin.

We add a disk by specifying a center and radius (and optionally a name):
```@example tut1
d = Disk([0,0], 0.8)
push!(bd, d)
# Make the structure required:
billiard = Billiard(bd)
```
To make sure the billiard looks as you would expect, use the function `plot_billiard(bd)`. Create a particle inside that billiard and evolve it:
```@example tut1
using PyPlot
plot_billiard(billiard)
ω = 0.5
p = randominside(billiard, ω)
xt, yt, vxt, vyt, t = construct(evolve!(p, billiard, 100)...)
plot(xt, yt)
plot_particle(p)
savefig("tut1.svg"); nothing # hide
```
![](tut1.svg)


The billiard table now works for straight or magnetic propagation.
To expand this to ray-splitting you have to use ray-splitting `Obstacle`s ([see the tutorial on Ray-Splitting](/tutorials/ray-splitting)).
Additional information on how to define your own `Obstacle` sub-type is given in the tutorial on [Defining your own Obstacles](/tutorials/own_obstacle).

If you make *any* billiard system that you think is cool and missing from this package, you are more than welcome to submit a PR extending the Standard Billiards Library with your contribution!

## Obstacle order
!!! info "Obstacle order."
    The order that the obstacles are given to the constructor is important for the
    function [`boundarymap`](@ref). For any other functionality it is irrelevant.

## Convex Billiards
These 2 types of walls used by `DynamicalBilliards` that behave differently during
evolution:
  1. [`InfiniteWall`](@ref) : This wall is not actually infinite. It has a starting and ending
     position. However, when the collision time is calculated, this wall is assumed
     to be a line (i.e. *infinite*). This is absolutely fine, as long as the
     billiards used are [*convex* polygons](https://en.wikipedia.org/wiki/Convex_polygon).
     [`SplitterWall`](@ref), [`PeriodicWall`](@ref) and [`RandomWall`](@ref) behave like `InfiniteWall` during evolution.
  2. `FiniteWall` : This wall is indeed finite in every sense of the word. This
     means that during collision time estimation, if the collision point that was
     calculated lies *outside* of the boundaries of the `FiniteWall`, then the
     returned collision time is `Inf` (no collision). `FiniteWall` is slower
     than `InfiniteWall` for that reason.


If you wish to create a billiard table that you know will be convex, you should
then use `InfiniteWall`s for faster evolution.
Notice that using [`escapetime`](@ref) requires
at least one `FiniteWall` with field `isdoor=true`.

## Obstacle Library
This is the list of `Obstacle`s you can use when creating your own billiard.

```@docs
Disk
RandomDisk
Antidot
Semicircle
InfiniteWall
RandomWall
PeriodicWall
SplitterWall
FiniteWall
```
---
In addition, `translate` is a helpful function:
```@docs
translate
```
