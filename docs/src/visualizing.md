All plotting functionality of `DynamicalBilliards` lies within a few well-defined functions that use the `PyPlot` package to plot aspects of the system on the current PyPlot figure.


**REWORK THIS: INSTALL PYPLOT AND BUILD PROPER DODCSTRINGS**

All plotting functions are brought into scope when `using PyPlot`. The functions are:
```@docs
plot_obstacle!
plot_particle!
plot_cyclotron!
plot_billiard
animate_evolution
```

## Examples

### Plotting some obstacles

For example:
```julia
using DynamicalBilliards, PyPlot
bd = billiard_sinai(0.3)
# Plot disk:
plot_obstacle!(bd[5])
# Plot left wall:
plot_obstacle!(bd[1])
# Plot right wall with different settings:
plot_obstacle!(bd[2]; linewidth = 3.0, linestyle = "dashed", color = (1.0, 0.5, 0.5))
# Set limits for the purpose of the tutorial
xlim(-0.1, 1.1); ylim(-0.1, 1.1)
```
will plot something like this:

![Visualizing tutorial 1](http://i.imgur.com/lrDStnP.png)

If you want to quickly plot the entire billiard with default parameters, simply use the function `plot_billiard(bd)`:
```julia
bd = billiard_polygon(6, 1)
a = Antidot([0.0,0.0], 0.5)
bd = Billiard(bd.obstacles..., a)
plot_billiard(bd)
```
which will plot something like this:

![Visualizing tutorial 2](http://i.imgur.com/46AomXm.png)

`plot_billiard()` also sets up the axis to have equal aspect ratio and sets up the axis limits to be just large enough to contain the entire billiard.



### Plotting particles

Following the above example, we create and plot a particle using the function `plot_particle`:
```julia
p = randominside(bd)
plot_particle(p)
# Plot one more particle with purple color,
# pentagon shape and bigger size (default is s=30):
p2 = randominside(bd)
plot_particle(p2; color=(0.5, 0, 0.8), marker="p", s=60.0)
```
which should give you something like this (notice that the particle position and direction are random):

![Visualizing tutorial 3](http://i.imgur.com/8a4ajfA.png)

## Color conventions
The default plotting settings have been chosen for maximum clarity and consistency. The color conventions followed are:
* Particles are black.
* Particle orbits/trajectories are blue.
* Reflecting obstacles (e.g. `Disk`, `FiniteWall` etc.) are green.
* Randomly reflecting obstacles (e.g. `RandomDisk` or `RandomWall`) are purple.
* Ray-splitting obstacles are red with dashed linestyle.
* Periodicity enforcing obstacles are yellow with dotted linestyle
  (if and when plotted).
* Doors (`InfiniteWall` with `isdoor=true`) are plotted with alternating black and
  cyan dashed lines.

## Animating the motion of a particle

The function `animate_evolution` is provided to animate the evolution of a particle from collision to collision:
```julia
animate_evolution(p, bd, colnumber[, ray-splitter]; kwargs...)
```

Arguments:
  * `p::AbstractParticle` : The particle to be evolved (gets mutated!).
  * `bd::Billiard` : The billiard.
  * `colnumber::Int` : Number of collisions to evolve the particle for.
  * `ray-splitter::Dict{Int, Any}` : (Optional) Ray-splitting dictionary
      that enables ray-splitting processes during evolution.

Keyword Arguments:
  * `newfig = true` : Creates a new figure at the function call, and plots
    the billiard in that figure.
  * `sleeptime` : Time passed to `sleep()` between each collision.
  * `col_to_plot` : How many previous collisions are shown during the animation.
  * `particle_kwargs` : Either a Dict{Symbol, Any} or a vector of Tuple{Symbol, Any}.
    Keywords passed into `plot_particle()`.
  * `orbit_kwargs` : Either a Dict{Symbol, Any} or a Vector of Tuple{Symbol, Any}.
    Keywords passed into `PyPlot.plot()` which plots the orbit of the particle
    (`line` object).
  * `savefigs::Bool` : If `true` save .png figures of each frame of the animation
    A direct movie (like creating a .mp4) of the animation cannot be made automatically,
    since the animation process mutates the particle.
  * `savename` : Name (*including path*) of the figures to be produced. The ending
    "\_i.png" will be attached to all figures.

The function returns `a, b, c`. Do `a[:remove](), b[:remove](), c[:remove]()` to clear
the particle out of the figure.

Automatic output into an animated image (e.g. ".gif" format) is not yet supported.

Let's animate a particle inside a simple pentagon with magnetic field:

```julia
bd = billiard_polygon(5, 1)
a = Disk([0.0,0.0], 0.4)
bd = (bd.obstacles..., a)
p = randominside(bd, 1.0)

savedir = tempdir()
animate_evolution(p, bd, 50; savefigs = true, savename = savedir)
```

This code produced 50 ".png" images which were later combined
into a single ".gif" animation:

![Visualizing Animation 1](http://i.imgur.com/UyiW2N2.gif)

## Periodic Billiards
In order to plot periodic billiards, you have need to call a different method of
[`plot_billiard`](/basic/library/#DynamicalBilliards.plot_billiard), since now you
also have to specify the limits of plotting. The
methods provided are:
```julia
plot_billiard(bd, xmin, ymin, xmax, ymax)
plot_billiard(bd, xt::Vector{T}, yt::Vector{T})
```
The last one conveniently plots the combo of particle-trajectory and
periodic-billiard taking care of all the details internally. Give the keyword
`plot_orbit = false` if you do not want to plot the orbit defined by `(xt, yt)`.

For example, the following code
```julia
using DynamicalBilliards, PyPlot
r = 0.25
bd = billiard_rectangle(2, 1; setting = "periodic")
d = Disk([0.5, 0.5], r)
d2 = Disk([1.5, 0.5], r/2)
bd = Billiard(bd.obstacles..., d, d2)
p = randominside(bd)
xt, yt, vxt, vyt, t = construct(evolve!(p, bd, 50)...)
plot_billiard(bd, xt, yt)
plot_particle(p)
```
will produce something like this:
![Periodic Billiard plot](http://i.imgur.com/rOpU7sl.png)

## Boundary Map plots
```@docs
plot_boundarymap
```
