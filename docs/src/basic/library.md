Below you find all the docstrings of all exported names of `DynamicalBilliards.jl` in convenient groups.

## Standard Billiards

```@docs
billiard_rectangle
billiard_sinai
billiard_polygon
```


## Particles

```@docs
Particle
MagneticParticle
magnetic2standard
standard2magnetic
cyclotron
```

## Obstacles

```@docs
Obstacle
Disk
Antidot
FiniteWall
PeriodicWall
SplitterWall
normalvec
distance
randominside
```

## Propagation

```@docs
resolvecollision!
collisiontime
propagate!
evolve!
construct
isphysical
```

## Ray-splitting

```@docs
isphysical
acceptable_raysplitter
supports_raysplitting
```

## Visualization
```@docs
plot_obstacle
plot_billiard
plot_particle
plot_cyclotron
animate_evolution
```

