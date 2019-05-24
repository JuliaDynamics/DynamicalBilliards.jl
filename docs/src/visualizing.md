All plotting functionality of `DynamicalBilliards` lies within a few well-defined functions that use the `PyPlot` package to plot aspects of the system.

All plotting functions are brought into scope when `using PyPlot`. The functions are:
```@docs
DynamicalBilliards.plot(::Obstacle)
DynamicalBilliards.plot(::Billiard)
DynamicalBilliards.plot(::AbstractParticle, ::Bool)
animate_evolution
```

## Examples

### Plotting Obstacles with keywords
```@example obstacles
using DynamicalBilliards, PyPlot

bd = billiard_sinai()

figure()
plot(bd[2]);
plot(bd[4], color = "blue", linestyle = "dotted", lw = 5.0);
plot(bd[1], facecolor = "yellow", edgecolor = "black");
savefig("rand_obstacles.svg"); nothing # hide
```
![](rand_obstacles.svg)

### Plotting a Billiard

```@example 8
using DynamicalBilliards
b = billiard_polygon(6, 1)
a = Antidot([0.0,0.0], 0.5)
bd = Billiard(b.obstacles..., a)
```

If you want to quickly plot the entire billiard with default parameters, simply use the function `plot(bd)`:

```@example 8
using PyPlot
plot(bd)
savefig("billiard_example.svg"); nothing # hide
```
![](billiard_example.svg)

`plot()` also sets up the axis to have equal aspect ratio and sets up the axis limits to be just large enough to contain the entire billiard.



### Plotting particles

Following the above example, we create and plot a particle using the function `plot`:
```@example 8
p = randominside(bd)
plot(p)
# Plot one more particle with purple color,
# pentagon shape and bigger size (default is s=30):
p2 = randominside(bd)
plot(p2; color=(0.5, 0, 0.8), marker="p", s=60.0)
p3 = randominside(bd, 2.0)
plot(p3, true; color=(0, 0, 0.8), marker="o", s=60.0)
savefig("particles_example.svg"); nothing # hide
```
![](particles_example.svg)

(notice that the particle position and direction are random)

## Color conventions
The default plotting settings have been chosen for maximum clarity and consistency. The color conventions followed are:
* Particles are black.
* Particle orbits use matplotlib's color cycle (first one is blue).
* Reflecting obstacles (e.g. `Disk`, `FiniteWall` etc.) are green.
* Randomly reflecting obstacles (e.g. `RandomDisk` or `RandomWall`) are purple.
* Ray-splitting obstacles are red with dashed linestyle.
* Periodicity enforcing obstacles are yellow with dotted linestyle
  (if and when plotted).
* Doors (`FiniteWall` with `isdoor=true`) are plotted with alternating black and cyan dashed lines.

## Animating the motion of a particle

The function [`animate_evolution`](@ref) is provided to animate the evolution of a particle from collision to collision.
Let's animate a particle inside a simple pentagon with magnetic field:

```julia
bd = billiard_polygon(5, 1)
a = Disk([0.0,0.0], 0.4)
bd = Billiard(bd.obstacles..., a)
p = MagneticParticle(0, -0.5, 0, 1.0)

animate_evolution(p, bd, 10.0; savename = "penta")
```

<video width="100%" height="auto" controls>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/billiards/penta.mp4?raw=true" type="video/mp4">
</video>

The cool thing about [`animate_evolution`](@ref) is that it can animate multiple
particles simultaneously! For example, here is a beautiful demonstration of chaos
in non-dispersive billiards:
```julia
bd = billiard_stadium()
N = 20
cs = [(i/N, 0, 1 - i/N, 0.5) for i in 1:N]
ps = [Particle(1, 0.6 + 0.0005*i, 0) for i in 1:N]

animate_evolution(ps, bd, 10.0; colors = cs, tailtime = 1.5, savename = "disperse")
```

<video width="100%" height="auto" controls>
<source src="https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/billiards/disperse.mp4?raw=true" type="video/mp4">
</video>




## Periodic Billiards
In order to plot periodic billiards, you have need to call a different method of
[`plot`](@ref), since now you
also have to specify the limits of plotting. The
methods provided are:
```julia
plot(bd, xmin, ymin, xmax, ymax)
plot(bd, xt::Vector{T}, yt::Vector{T})
```
The last one conveniently plots the combo of particle-trajectory and
periodic-billiard taking care of all the details internally. Give the keyword
`plot_orbit = false` if you do not want to plot the orbit defined by `(xt, yt)`.

For example, the following code
```@example 8
using DynamicalBilliards, PyPlot
r = 0.25
bd = billiard_rectangle(2, 1; setting = "periodic")
d = Disk([0.5, 0.5], r)
d2 = Ellipse([1.5, 0.5], r, 2r/3)
bd = Billiard(bd.obstacles..., d, d2)
p = Particle(1.0, 0.5, 0.2)
xt, yt, vxt, vyt, t = timeseries!(p, bd, 10)
plot(bd, xt, yt)
plot(p)
savefig("rectperiodic.svg"); nothing # hide
```
![](rectperiodic.svg)

And, you can also periodically plot billiards with hexagonal periodicity. Only give
the keyword argument `hexagonal = true` to [`plot`](@ref). As an example:

```@example 8
bd = billiard_hexagonal_sinai(0.3, 1.50; setting = "periodic")
d = Disk([0.7, 0], 0.2)
d2 = Antidot([0.7/2, 0.65], 0.35)
bd = Billiard(bd..., d, d2)

p = MagneticParticle(-0.5, 0.5, π/5, 1.0)

xt, yt = timeseries(p, bd, 10)

plot(bd, xt, yt; hexagonal = true)
savefig("hexperiodic.svg"); nothing # hide
```
![](hexperiodic.svg)


## Boundary Map plots
```@docs
plot_boundarymap
```
Examples are shown in the [phase spaces](basic/phasespaces.md) section.
