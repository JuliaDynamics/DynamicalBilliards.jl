Below you find all the docstrings of all exported names of `DynamicalBilliards.jl` in convenient groups.

## Standard Billiards

```@docs
billiard_rectangle
billiard_sinai
billiard_polygon
billiard_julia
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
Circular
Disk
RandomDisk
Antidot
Wall
FiniteWall
RandomWall
PeriodicWall
SplitterWall
normalvec
distance
randominside
```

## Propagation

```@docs
resolvecollision!
relocate!
specular!
lct
periodicity!
collisiontime
propagate!
evolve!
construct
realangle
```

## Ray-splitting

```@docs
isphysical
acceptable_raysplitter
supports_raysplitting
```

## Visualization
All visualization-related functions are contained in the `DynamicalBilliardsPlotting` package.
Here are their documentation strings:

```julia
plot_obstacle(obst::Obstacle; kwargs...)
```
Plot given obstacle on the current `PyPlot` figure. The default arguments for
each type of obstacle have been chosen for maximum clarity and
consistency. The `kwargs...` given by the user are keywords passed directly into PyPlot's
constructors. For `Wall` obstacles, kwargs are passed into `PyPlot.plot()`. For
`Disk` obstacles, kwargs are passed into `PyPlot.plt[:Circle]()`.

```julia
plot_particle(p::AbstractParticle; use_cell=true, kwargs...)
```
Plot given particle on the current `PyPlot` figure. Optionally use `p.current_cell` for
the particle's position. Given `kwargs...` are passed onto `PyPlot.scatter()`.
The particle is represented as a small ball (`PyPlot.scatter()`) and a small arrow (`PyPlot.quiver()`).
All `kwargs...` are given to `scatter()` but if a keyword argument `color` is given,
it is also passed to `quiver()`.

```julia
plot_cyclotron(p::MagneticParticle; use_cell=true, kwargs...)
```
Plot the circle traced by the free particle motion. Optionally use `p.current_cell` for
the particle's position. The user provided `kwargs...` are passed onto `PyPlot.plt[:Circle]()`.

```julia
plot_billiard(bt::Vector{Obstacle})
```
Plot all obstacles in `bt` using the default arguments, set
`xlim` and `ylim` to be 10% larger than `cellsize` and
set the axis aspect ratio to equal.

```julia
plot_billiard(bt, xmin, ymin, xmax, ymax)
```
Plot the given (periodic) billiard `bt` on the current PyPlot figure, repeatedly
plotting from `(xmin, ymin)` to `(xmax, ymax)`. Only works for rectangular billiards.

```julia
plot_billiard(bt, xt::Vector{Float64}, yt::Vector{Float64}; plot_orbit = true)
```
Plot the given (periodic) billiard `bt` along with a particle trajectory defined
by `xt` and `yt`, on the current PyPlot figure. Only works for rectangular billiards.
Sets limits automatically. Set the keyword argument `plot_orbit = false` to not
plot the orbit defined by `(xt, yt)`.

```julia
billiard_julia(; plotit = true)
```
Return the awesome "Julia-logo" billiard shown in the introduction
of DynamicalBilliards.jl. By default it also plots the billiard in a new
`PyPlot.figure()` using the correct colors.

### Animations

```julia
animate_evolution(p, bt, colnumber[, ray-splitter]; kwargs)
```
Animate the evolution of the particle, plotting the orbit from collision to collision.

Notice the difference with `evolve!()`: No time is given here; instead a number of
collisions is passed.

#### Arguments
* `p::AbstractParticle` : Either standard or magnetic.
* `bt::Vector{Obstacle}` : The billiard table.
* `colnumber::Int` : Number of collisions to evolve the particle for.
* `ray-splitter::Dict{Int, Any}` : (Optional) Ray-splitting dictionary
  that enables ray-splitting processes during evolution.
#### Keyword Arguments
* `sleeptime` : Time passed to `sleep()` between each collision.
* `col_to_plot` : How many previous collisions are shown during the animation.
* `savefigs` : Save .png figures to enable the creation of animation afterwards.
  **WARNING:** currently the .gif production has to be made by the user!
* `savename` : Name (**including path!**) of the figures to be produced. The ending
  "\_i.png" will be attached to all figures.
* `particle_kwargs` : Either a Dict{Symbol, Any} or a vector of Tuple{Symbol, Any}.
  Keywords passed into `plot_particle()`.
* `orbit_kwargs` : Either a Dict{Symbol, Any} or a Vector of Tuple{Symbol, Any}.
  Keywords passed into `PyPlot.plot()` which plots the orbit of the particle.
