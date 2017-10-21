All plotting functionality of `DynamicalBilliards` lies within a few well-defined functions that use the `PyPlot` package to plot aspects of the system on the current PyPlot figure.

*Remember to use* `DynamicalBilliards.enableplotting()` *to bring the plotting functions into scope!*

The functions are:
```julia
plot_obstacle(obst::Obstacle; kwargs...)
plot_particle(p::AbstractParticle; use_cell=true, kwargs...)
plot_cyclotron(p::MagneticParticle; use_cell=true, kwargs...)
plot_billiard(bt::Vector{Obstacle})
plot_billiard(bt, xt::Vector, yt::Vector; plot_orbit = true)
animate_evolution!(p, bt, colnumber[, ray-splitter]; kwargs)
billiard_julia(; plotit = true)
```
and each has a detailed documentation string.

## Plotting the Billiard Table

The function `plot_obstacle(obst::Obstacle; kwargs...)` plots the given obstacle on the current PyPlot figure. The `kwargs...` are keywords passed directly into `PyPlot`'s constructors (like e.g. `linewidth = 2.0`).
For example:
```julia
using DynamicalBilliards, PyPlot
DynamicalBilliards.enableplotting()
bt = billiard_sinai(0.3)
# Plot disk:
plot_obstacle(bt[5])
# Plot left wall:
plot_obstacle(bt[1])
# Plot right wall with different settings:
plot_obstacle(bt[2]; linewidth = 3.0, linestyle = "dashed", color = (1.0, 0.5, 0.5))
# Set limits for the purpose of the tutorial
xlim(-0.1, 1.1); ylim(-0.1, 1.1)
```
will plot something like this:

![Visualizing tutorial 1](http://i.imgur.com/lrDStnP.png)

If you want to quickly plot the entire billiard table without changing the settings, simply use the function  `plot_billiard(bt)`:
```julia
bt = billiard_polygon(6, 1)
a = Antidot([0.0,0.0], 0.5)
push!(bt, a)
plot_billiard(bt)
```
which will plot something like this:

![Visualizing tutorial 2](http://i.imgur.com/46AomXm.png)

`plot_billiard()` also sets up the axis to have equal aspect ratio and sets up the axis limits to be just large enough to contain the entire billiard.



## Plotting particles

Following the above example, we create and plot a particle using the function `plot_particle`:
```julia
p = randominside(bt)
plot_particle(p)
# Plot one more particle with purple color,
# pentagon shape and bigger size (default is s=30):
p2 = randominside(bt)
plot_particle(p2; color=(0.5, 0, 0.8), marker="p", s=60.0)
```
which should give you something like this (notice that the particle position and direction are random):

![Visualizing tutorial 3](http://i.imgur.com/8a4ajfA.png)

## Color conventions
The default plotting settings have been chosen for maximum clarity and consistency. The color conventions followed are:
* Particles are black.
* Particle orbits/trajectories are blue.
* Reflecting obstacles (e.g. `Disk`, `FiniteWall` or `InfiniteWall`) are green.
* Randomly reflecting obstacles (e.g. `RandomDisk` or `RandomWall`) are yellow.
* Ray-splitting obstacles are red with dashed linestyle.
* Periodicity enforcing obstacles are purple with dotted linestyle
  (if and when plotted).
* Doors (`InfiniteWall` with `isdoor=true`) are plotted with alternating black and
  cyan dashed lines.

## Animating the motion of a particle

The function `animate_evolution!` is provided to animate the evolution of a particle from collision to collision:

```julia
animate_evolution!(p, bt, colnumber[, ray-splitter]; kwargs...)
```
Animate the evolution of the particle, plotting the orbit from collision to collision.

### Arguments
  * `p::AbstractParticle` : The particle to be evolved (gets mutated!).
  * `bt::Vector{Obstacle}` : The billiard table.
  * `colnumber::Int` : Number of collisions to evolve the particle for.
  * `ray-splitter::Dict{Int, Any}` : (Optional) Ray-splitting dictionary
      that enables ray-splitting processes during evolution.
### Keyword Arguments
  * `newfig = true` : Creates a new figure at the function call, and plots
    the billiard table in that figure.
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

Automatic output into an animated image (e.g. ".gif" format) is not possible, because
the function mutates its argument. However, `animate_evolution` gives users the possibility
to save each produced figure in order to merge as an animation using an external tool.

Let's animate a particle inside a simple pentagon with magnetic field:

```julia
bt = billiard_polygon(5, 1)
a = Disk([0.0,0.0], 0.4)
push!(bt, a)
p = randominside(bt, 1.0)

savedir = tempdir()
animate_evolution!(p, bt, 50; savefigs = true, savename = savedir)
```

This code produced 50 ".png" images which were later mixed (using e.g. [gifmaker](www.gifmaker.me)) into a single ".gif" animation:

![Visualizing Animation 1](http://i.imgur.com/UyiW2N2.gif)

## Periodic Billiards
In order to plot periodic billiards, you have need to call a different method of
[`plot_billiard`](/basic/library/#DynamicalBilliards.plot_billiard), since now you
also have to specify the limits of plotting. The
methods provided are:
```julia
plot_billiard(bt, xmin, ymin, xmax, ymax)
plot_billiard(bt, xt::Vector{T}, yt::Vector{T})
```
The last one conveniently plots the combo of particle-trajectory and
periodic-billiard taking care of all the details internally. Give the keyword
`plot_orbit = false` if you do not want to plot the orbit defined by `(xt, yt)`.

For example, the following code
```julia
using DynamicalBilliards
DynamicalBilliards.enableplotting()
r = 0.25
bt = billiard_rectangle(2, 1; setting = "periodic")
d = Disk([0.5, 0.5], r)
d2 = Disk([1.5, 0.5], r/2)
push!(bt, d, d2)
p = randominside(bt)
xt, yt, vxt, vyt, t = construct(evolve!(p, bt, 50)...)
plot_billiard(bt, xt, yt)
plot_particle(p)
```
will produce something like this:
![Periodic Billiard plot](http://i.imgur.com/rOpU7sl.png)

Animations for periodic billiards are not supported yet.
