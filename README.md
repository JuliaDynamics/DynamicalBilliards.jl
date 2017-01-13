# DynamicalBilliards.jl
A Julia package for dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction. ![Example animation](https://github.com/Datseris/DynamicalBilliards.jl/blob/master/images/plot_example.gif "Evolution of particle in a magnetic field.")

*(the code that generated this animation is shown at the end of this README)*

The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:
* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard. 
* **Ray-splitting Billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiard systems that have already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : functions for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more. Also includes animation related content.

NOTICE: This package does not support collision between particles. All particles are considered point-particles for all simulations offered by `DynamicalBilliards.jl`.

Documentation: coming soon.

*(all exported names of DynamicalBilliards.jl have detailed docstrings. Use `?` when in doubt)*

---
### Installation
This package is not yet registered. Use `Pkg.clone("https://github.com/Datseris/DynamicalBilliards.jl")` in order to install it.
## Basic Usage
The simplest usage of this package revolves around a single function: 
```julia
evolve!(p::Particle, bt::Vector{Obstacle}, total_time)
```
The function evolves a particle `p` inside a billiard table `bt` for a given amount of time `total_time`, while taking care of all the details internally. 

The first step is to define the billiard table `bt`, which is the system the particle `p` will propagate in. A billiard table is simply a collection (`Vector`) of `Obstacle`s. The most convenient way is to use one of the pre-defined billiard tables offered by the package. For example, let's create a periodic Sinai billiard with disk radius of 0.3 and with one side of length 2 and one of length 1:
```julia
bt = billiard_sinai_periodic(0.3, 2.0, 1.0)                                                    
```
(for more information about defining billiard tables see the official documentation here)

Afterwards, you would want to create a particle inside that billiard system. For that, the function `randominside(bt::Vector{Obstacle})` is provided. This function returns a particle with random initial conditions inside the billiard table, while making sure that it is always in the allowed region of the billiard table.
```julia
p = randominside(bt)
```
Now you are ready to evolve this particle:
```julia
t, poss, vels = evolve!(p, bt, 1000.0)
```
The return values of the `evolve!()` function need some brief explaining: As noted by the "!" at the end of the function, it changes its argument `p` (specifically, it updates almost all the fields of `p`).
Most importantly however, this function also returns the main output expected by a billiard
system. This output is a tuple of three vectors:
* `t::Vector` : Collision times.
* `poss::Vector{SVector{2}}` : Positions during collisions.
* `vels:: Vector{SVector{2}})` : Velocities **exactly after** the collisions.
The time `t[i]` is the time necessary to reach state `poss[i], vels[i]` starting from the
state `poss[i-1], vels[i-1]`. That is why `t[1]` is always 0 since `poss[0], vels[0]` are
the initial conditions.

If this output is not convenient for you, the function `construct(t, poss, vels, dt=0.1*one(T))` is provided, which constructs the (continuous) timeseries of the position and velocity, as well as the time-vector, when given the main output of `evolve!()`:
```julia
xt, yt, vxt, vyt, ts = construct(t, poss, vels)
```
or, if you want to use the fancy ellipsis operator, you can do:
```julia
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0)...)
```

## Magnetic Propagation
The are only two differences between magnetic and straight propagation. Firstly, the particle type is not `Particle` anymore, but `MagneticParticle`. The latter has an extra field called `omega` which is the cyclic frequency of rotation (equivalently, the angular velocity). In order to create a `MagneticParticle` (without using the constructors), you simply provide this extra argument to the `randominside()` function:
```julia
ω = 0.5
p = randominside(ω, bt)
typeof(p) # MagneticParticle
p.omega   # 0.5
```
To propagate the particle you use the same functions
```
t, poss, vels = evolve!(p, bt, 1000.0)
xt, yt, vxt, vyt, ts = construct(ω, t, poss, vels)
# or equivalently: 
xt, yt, vxt, vyt, ts = construct(ω, evolve!(p, bt, 1000.0)...)
```
As you can see, the second difference is that the additional argument of the angular velocity must also be provided to the `construct()` function, in order for it to construct circular motion instead of straight motion between collisions. (Note: the `ω` argument is always given as the first argument, for consistency)

## Visualizing
*(all plotting in* `DynamicalBilliards` *is currently done through the* `PyPlot` *package. In a future update a switch will happen towards* `Plots.jl` *)*

The functions `plot_obstacle(o::Obstacle; kwargs...)`, `plot_billiard(bt::Vector{Obstacle}; kwargs...)` and `plot_particle(p::AbstractParticle; kwargs...)` are provided in order to plot the respective elements **on the current PyPlot figure**. In order to animate the evolution of a particle in a billiard table, use the function:
```julia
plot_evolution(p::AbstractParticle, bt::Vector{Obstacle}, colnumber = 50;
               sleeptime = 0.5, col_to_plot = 5, color = (0,0,1), savefigs = false, savename = "")
```
which propagates the particle from collision to collision up to a total of `colnumber` collisions. Then, it draws the collisions, drawing always only the last `col_to_plot` collisions with orbit color `color`.  `sleeptime` signals the waiting time between each plot update.

Direct animation saving is not supported yet. However, optionally, you could save each figure of the animation using `savefigs = true` and `savename = "path/to/your/folder/figurename"`. A total of `colnumber` figures will be created, ending with `_#.png`.

Be sure to first call `plot_billiard` before calling `plot_evolution`.

The example .gif shown in the beginning of this README, was generated simply with the code:
```julia
using DynamicalBilliards
using PyPlot

bt = billiard_rectangle(1.5, 1.0)
d1 = Disk([0.45, 0.6], 0.3, "Upper-left Disk")
d2 = Disk([1.1, 0.3], 0.15, "Lower-right Disk")
d3 = Disk([1.2, 0.8], 0.1, "Small Disk")
w1 = FiniteWall([0.0, 0.4], [0.6,0.0], [0.4,0.6], "Diagonal")
push!(bt, d1, d2, d3, w1)
ω = 2.0
p = randominside(ω, bt)

plot_billiard(bt)
axis("off")
tight_layout()
xlim(0,1.5)
ylim(0,1.0)

sname = "C:\\***\\example"
plot_evolution(p, bt, 200;
sleeptime = 0.1, col_to_plot = 4, savefigs = true, savename = sname)
```
Afterwards the outputed .png files where merged into a single .gif externally.

A full overview of all plotting procedures offered by this package is coming soon.
