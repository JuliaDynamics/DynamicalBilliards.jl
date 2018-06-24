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
Boundary maps can be obdained with the high level function
```@docs
boundarymap
```
---
For example, take a look at boundary maps of the mushroom billiard, which is known to have a mixed phase space:
```@example coords
using DynamicalBilliards

bd = billiard_mushroom()

n = 100 # how many particles to create

ξς, φς, ις = boundarymap(bd, 10000, n)

using PyPlot # enables plot_boundarymap function

colors = ["C$(rand(1:9))" for i in 1:n] # random colors

figure()
plot_boundarymap(ξς, φς, ις, color = colors)
savefig("boundarymap.svg"); nothing # hide
```
![](boundarymap.svg)

<!-- ![Boundary map](https://i.imgur.com/RO9UZa9.png) -->

And of course similarly for magnetic fields
```@example coords
ξς, φς, ις = boundarymap(bd, 10000, n, 1.0) # angular velocity last argument
figure()
plot_boundarymap(ξς, φς, ις, color = colors)
savefig("boundarymapmag.svg"); nothing # hide
```
![](boundarymapmag.svg)

<!-- ![Boundary map with magnetic field](https://i.imgur.com/YoW1FVD.png) -->

## Phase Space Portions
It is possible to compute the portion of phase space covered by a particle as it
is evolved in time. We have two methods, one for the "boundary" coordinates (2D space)
and one for the "real" coordinates (3D space):
```@docs
boundarymap_portion
phasespace_portion
```
For example, for mushroom billiards the ratio of the chaotic-to-total phase space is known **analytically** for both the full 3D [1] space as well as the boundary 2D [2] space:
$$
formulas from Lukas Thesis
$$
We can easily confirm those formulas:
```@example phasespace
using DynamicalBilliards

t = 1000000.0
l = 1.0; r = 1.0; w = 0.4

bd = billiard_mushroom(l, w, r)

p = MushroomTools.randomchaotic(l, w, r)

ratio, dic = boundarymap_portion(bd, t, p, 0.01)
trueratio = MushroomTools.g_c_2D(l,w,r)
println("2D numeric - theory: $(abs(ratio - trueratio))")

ratio = phasespace_portion(bd, t, p, 0.01)
trueratio = MushroomTools.g_c_3D(l,w,r)
println("3D numeric - theory: $(abs(ratio - trueratio))")
```
Of course, increasing evolution time or decreasing boxsize will bring higher accuracy.

## Chaotic vs. Regular boundary map animation
A simple 1-2 sentences summary
```
script that produces it
```

animation file.

To be done by Lukas Hupe.

---

## References

[1] : [A. H. Barnett & T. Betcke, *Quantum mushroom billiards*, Chaos, 17(4) (20017)](https://doi.org/10.1063/1.2816946)

[2] : Lukas Hupe, B.Sc. Thesis (2018), *to be published*
