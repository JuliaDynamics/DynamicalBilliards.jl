# Setting up your Billiard Table
A billiard table `bt` is a vector of Obstacles: `bt::Vector{Obstacle}`. 
The abstract Type `Obstacle` is the supertype of all objects that a particle may collide with.

There are some premade functions that construct well-known billiards, like the periodic Sinai billiard.
You can find all of them at the Standard Billiards page.

To create a custom billiard, you start with an empty Vector:
```
bt = Obstacle[]
```
and then you create your obstacles one by one and add them to it. All obstacles that are already defined in the package
can be found at the Obstacles page. For the example of this page, we will create a hexagonal billiard with a disk in the middle.

The first step is to define the six walls of the billiard table. A `FiniteWall` object needs to be supplemented with a start point, an end point, a normal vector and, optionally, a name. The vertex points of a regular hexagon of radius `r` are given by the formula:
```math
(x,y) = \left( r\cos\left(\frac{2\pi i}{6}\right), r\cos\left(\frac{2\pi i}{6}\right) \right) \quad \text{for i $\in$ \{1,...,6\}}
```
To create each wall object, we will implement the following loop, choosing a size of 2.0:
```julia
hexagon_vertex = (r) -> [ [r*cos(2π*i/6), r*sin(2π*i/6)] for i in 1:6]
hexver = hexagon_vertex(2.0)

for i in eachindex(hexver)
  starting = hexver[i]
  ending = hexver[mod1(i+1, length(hexver))]
  w = ending - starting
  normal = [-w[2], w[1]]
  wall = FiniteWall(starting, ending, normal, "wall $i")
  push!(bt, wall)
end
```
The `normal` vector of a `Wall` obstacle is necessary to be supplemented by the user because it must point towards where the particle is expected to come from. If `w` is the vector (wall) pointing from start- to end-point then the vector `[-w[2], w[1]]` is pointing to the left of `w` and the vector `[w[2], -[w1]]` is pointing to the right. Both are normal to `w`, but you have to know which one to pick. In this case this is very easy, since the normal has to simply point towards the origin.

We add a disk by specifying a center and radius (and optionally a name): 
```julia
d = Disk([0,0], 0.8)
push!(bt, d)
```
To make sure the billiard looks as you would expect, use the function `plot_billiard(bt)`. Create a particle inside that billiard and evolve it while also plotting it:
```julia
plot_billiard(bt)
ω = 0.5
p = randominside(p, ω)
plot_evolution(p, bt, 100)
```
    
The billiard table now works for straight or magnetic propagation. To expand this to ray-splitting you have to use ray-splitting Obstacles (see the tutorial on Ray-Splitting). Additional information on how to define your own `Obstacle` sub-type is given in the tutorial on Defining your own Obstacles.

If you make any billiard system that you think is common and missing from this package, you are more than welcome to submit a PR extending the `StandardBilliards.jl` library with your contribution!
    
