# Internals

This page is **not** part of the public API defined by `DynamicalBilliards`. Consider it something like a *developer's guide*.

## Implementation
Before talking about the low level methods that enable everything to work nicely together, let's talk about how this package works.

Firstly one defines a [`Billiard`](@ref) and optionally some [`RaySplitter`](@ref) instances. Then one creates a particle inside the defined billiard. The algorithm for the propagation of a particle is the following:

1. Calculate the [`collision`](@ref) of the particle with **all** obstacles in the billiard.
2. Find the collision that happens first (in time), and the obstacle corresponding to that.
3. [`DynamicalBilliards.relocate!`](@ref) the particle, and ensure that it is **inside** the billiard. This means that [`DynamicalBilliards.distance`](@ref) between particle and obstacle is either positive or close to machine precision.
4. (Optionally) check if there is transmission for ray-splitting: `T(φ) > rand()`
  * If yes, perform the ray-splitting algorithm (see the [ray-splitting](ray-splitting) page).
  * If not, then [`DynamicalBilliards.resolvecollision!`](@ref) of the particle with the obstacle (specular or periodic conditions).

5. Continue this loop for a given amount of time.

Notice that the [`DynamicalBilliards.relocate!`](@ref) step is *very* important because it takes care that all particles remain inside the billiard.

The exposed [`bounce!`](@ref) function bundles steps 1-4 together.

### Where is "inside"?
If for some reason (finite numeric precision) a particle goes outside a billiard,
then it will escape to infinity. But what *is* inside?

"Inside" is defined on obstacle level by the function `distance`:
```@docs
DynamicalBilliards.distance
```
Notice that for very small negative values of distance, [`collision`](@ref) takes care of finite precision issues and does not return wrong collisions.


## Numerical Precision

All core types of `DynamicalBilliards` are parametrically constructed, with
parameter `T <: AbstractFloat`. This means that the fields of all particles and obstacles
contain numbers strictly of type `T`. You will understand why this choice happened
as you continue reading this paragraph.

The main concerns during evolution in a billiard table are:

1. The particle must never leak out of the billiard table. This is simply translated
   to the `distance` function being positive after any collision _and_ that `collision` takes care of extreme cases with very small (but negative) distance.
2. The collision time is never infinite, besides the cases of
   [pinned particles](physics/#pinned-particles) in a magnetic billiard.

These are solved with two ways:
1. After the next collision is computed, `relocate!` brings the particle to that point and calculates the `distance` with the colliding obstacle. If it is negative, it translates the particle's position by this distance, _along the normal vector_.
2. `collision` takes care of cases where the distance between particle and obstacle is less than `accuracy(::T)`. (This is necessary only for magnetic propagation, as for straight propagation checking the velocity direction with respect to the normal is always enough).

Adjusting the global precision of `DynamicalBilliards` is easy and can be done by choosing the floating precision you would like.
This is done by initializing your billiard table with parametric type `T`, e.g. `bd = billiard_sinai(Float16(0.3))`.
This choice will propagate to the entire `bd`, all particles resulting from [`randominside`](@ref), **as well as the entire evolution process**.

!!! danger "BigFloats"
    Evolution with `BigFloat` in `DynamicalBilliards` is on average
    3 to 4 orders of magnitude slower than with `Float64`.

---

## Collision Times
```@docs
collision
next_collision
```

## Non-Exported Internals

### Obstacle related

```@docs
DynamicalBilliards.normalvec
DynamicalBilliards.cellsize
```

### Propagation
```@docs
DynamicalBilliards.propagate!
DynamicalBilliards.resolvecollision!
DynamicalBilliards.relocate!
DynamicalBilliards.specular!
DynamicalBilliards.periodicity!
```

!!! warning "Cyclotron center is a field of `MagneticParticle`"
    For almost all operations involving a `MagneticParticle`, the center of
    the cyclotron is required. In order to compute this center only when it
    physically changes, we have made it a field of the `struct`.

    This means that after changing the position or velocity of the particle,
    this center must be changed by doing `mp.center = DynamicalBilliards.find_cyclotron(mp)`.
    The [`bounce!`](@ref) function takes care of that in the most opportune moment, but if you want to write your own specific low level function, do not forget this point!
