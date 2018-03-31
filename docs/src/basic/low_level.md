# Low Level API

## Implementation
Before talking about the low level methods that enable everything to work nicely
together, let's talk about how this package works.

Firstly one defines a [`Billiard`](@ref) and (if desired) the [ray-splitting dictionary](/tutorials/ray-splitting/#ray-splitter-dictionary). Then one creates a particle inside the defined billiard table. The algorithm for the propagation of a particle is the following:

1. Calculate the [`collisiontime`](@ref) of the particle with **all** obstacles in the billiard.
2. Find the smallest time, and the obstacle corresponding to that.
3. [`relocate!`](@ref) the particle, so that it is on the correct side of the obstacle to-be-collided with.
3. (Optionally) check if there is transmission for ray-splitting: `T(φ) > rand()`
    * If yes, perform the ray-splitting algorithm (not discussed here).
    * If not, then [`resolvecollision!`](@ref) of the particle with the obstacle.

4. Continue this loop for a given amount of time.

Notice that the [`relocate!`](@ref) step is *very* important because it takes care that all particles remain inside the billiard.

### Where is "inside"?
If for some reason (finite numeric precision) a particle goes outside a billiard,
then it will escape to infinity. But what *is* inside?

"Inside" is defined on obstacle level by the function `distance`:
```@docs
distance
```

---

## It's all about bounce!
The algorithm steps 1-3 described above are bundled in the following well-behaving function:
```@docs
bounce!
```
---
`bounce!` is the function used internally by all high-level functions, like [`evolve!`](@ref), [`poincaresection`](@ref), [`escapetime`](@ref), etc.

This is the function a user should use if they want to calculate other things besides what is already available in the high level API.


## Obstacle Library
```@docs
Obstacle
Circular
Disk
RandomDisk
Antidot
Semicircle
Wall
InfiniteWall
RandomWall
PeriodicWall
SplitterWall
FiniteWall
```
---
### Obstacle-related functions
```@docs
normalvec
cellsize
arclength
totallength
```

## Collision Times
```@docs
collisiontime
next_collision
realangle
```

## Propagation functions
```@docs
propagate!
resolvecollision!
relocate!
specular!
periodicity!
```


## Ray-splitting
```@docs
isphysical
acceptable_raysplitter
reset_billiard!
```
