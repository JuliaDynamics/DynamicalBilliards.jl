Below you find all the docstrings of all exported names of `DynamicalBilliards.jl` in convenient groups.



## Particles

```@docs
Particle
MagneticParticle
cyclotron
```

## Obstacles

```@docs
Obstacle
Circular
Disk
RandomDisk
Antidot
Wall
InfiniteWall
RandomWall
PeriodicWall
SplitterWall
normalvec
distance
randominside
```

## Standard Billiards

```@docs
billiard_rectangle
billiard_sinai
billiard_lorentz
billiard_polygon
billiard_hexagonal_sinai
billiard_raysplitting_showcase
```

## Propagation

```@docs
resolvecollision!
specular!
periodicity!
collisiontime
next_collision
propagate_pos
propagate!
evolve!
construct
```

## Ray-splitting

```@docs
isphysical
acceptable_raysplitter
reset_billiard!
```

## Visualization

Because the visualization related functions are loaded on-demand, use the Julia help mode
to read the documentation strings. With the command `DynamicalBilliards.enableplotting()`
the names that are brought into scope are:
```julia
plot_obstacle(obst::Obstacle; kwargs...)
plot_particle(p::AbstractParticle; use_cell=true, kwargs...)
plot_cyclotron(p::MagneticParticle; use_cell=true, kwargs...)
plot_billiard(bt::Vector{Obstacle})
plot_billiard(bt, xt::Vector, yt::Vector; plot_orbit = true)
billiard_julia(; plotit = true)
animate_evolution(p, bt, colnumber[, ray-splitter]; kwargs)
```
