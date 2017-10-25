Below you find all the docstrings of all exported names of `DynamicalBilliards.jl` in convenient groups.


## Propagation

```@docs
resolvecollision!
specular!
periodicity!
relocate!
collisiontime
next_collision
propagate_pos
propagate!
evolve!
construct
escapetime
```

## Particles

```@docs
Particle
MagneticParticle
cyclotron
```

## Obstacles

```@docs
Obstacle
Disk
RandomDisk
Antidot
Semicircle
InfiniteWall
RandomWall
PeriodicWall
SplitterWall
FiniteWall
normalvec
distance
randominside
```

## Standard Billiards

```@docs
billiard_rectangle
billiard_sinai
billiard_bunimovich
billiard_stadium
billiard_mushroom
billiard_lorentz
billiard_polygon
billiard_hexagonal_sinai
billiard_raysplitting_showcase
```

## Ray-splitting

```@docs
isphysical
acceptable_raysplitter
reset_billiard!
```
