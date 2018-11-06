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
```@example coords
using DynamicalBilliards, PyPlot

bd = billiard_mushroom()

n = 100 # how many particles to create
t = 200 # how long to evolve each one

bmap, arcs = parallelize(boundarymap, bd, t, n)

using PyPlot # enables plot_boundarymap function

colors = ["C$(rand(1:9))" for i in 1:n] # random colors

plot_boundarymap(bmap, arcs, color = colors)
tight_layout()
savefig("boundarymap.svg"); nothing # hide
```
![](boundarymap.svg)


And of course similarly for magnetic fields
```@example coords
bmap, arcs = parallelize(boundarymap, bd, t, n, 1.0)
plot_boundarymap(bmap, arcs, color = colors)
tight_layout()
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

---

For example, for mushroom billiards the ratio of the chaotic-to-total phase space is known **analytically** for both the full 3D [1] space as well as the boundary 2D [2] space:

```math
\begin{aligned}
  g_\text{c, 3D} &= \frac{2 \, \pi r^{2} - 4 \, r^{2} \arccos\left(\frac{w}{2 \, r}\right) + 4 \, l w + \sqrt{4 \, r^{2} - w^{2}} w}{2 \, {\left(\pi r^{2} + 2 \, l w\right)}}\\
  g_\text{c, 2D} &= \frac{\pi w + 2 \, w \arccos\left(\frac{w}{2 \, r}\right) + 4 \, l + 4 \, r - 2 \, \sqrt{4 \, r^{2} - w^{2}}}{2 \, {\left(\pi r + 2 \, l + 2 \, r\right)}}
\end{aligned}
```

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
Of course, increasing evolution time and decreasing boxsize will bring higher accuracy.

## References

[1] : A. H. Barnett & T. Betcke, *Quantum mushroom billiards*, [Chaos, 17(4) (20017)](https://doi.org/10.1063/1.2816946)

[2] : Lukas Hupe, B.Sc. Thesis (2018), *to be published*
