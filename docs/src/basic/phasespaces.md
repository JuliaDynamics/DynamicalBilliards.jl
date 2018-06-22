# Phase Spaces

## Coordinate Systems
Any billiard has two coordinate systems:
1. The "real" coordinates, i.e. the coordinates that specify the full, three-dimensional phase space: $x, y, \phi$.
2. The "boundary" coordinates, also known as Birkhoff coordinates, which instead reduce the continuous billiard into a map, by only considering the collision points. These coordinates are only two: $\xi, \sin(\phi_n)$, with $\xi$ being the parametrization of the arc length and $\phi_n$ being the angle as measured from the normal vector.

With `DynamicalBilliards` it is very easy to switch between the two coordinate
systems, using:
```@docs
to_bcoords
from_bcoords
arcintervals
```

## Boundary Maps
Boundary maps can be obtained with the high level function
```@docs
boundarymap
```
---
For example, take a look at boundary maps of the mushroom billiard, which is known to have a mixed phase space:
```julia
using DynamicalBilliards

bt = billiard_mushroom()

n = 100 # how many particles to create

ξς, φς, ις = boundarymap(bt, 10000, n)

using PyPlot # enables plot_boundarymap function

colors = ["C$(rand(1:9))" for i in 1:n] # random colors

figure()
plot_boundarymap(ξς, φς, ις, color = colors)
```
![Boundary map](https://i.imgur.com/RO9UZa9.png)

And of course similarly for magnetic fields
```julia
ξς, φς, ις = boundarymap(bt, 10000, n, 1.0) # angular velocity last argument
figure()
plot_boundarymap(ξς, φς, ις, color = colors)
```
![Boundary map with magnetic field](https://i.imgur.com/YoW1FVD.png)

## Phase Space Portions

```@docs
boundarymap_portion
phasespace_portion
```
