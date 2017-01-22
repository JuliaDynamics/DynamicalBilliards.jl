# Basic Usage
> *all names of* `DynamicalBilliards.jl` *have detailed docstrings. Use* `?` *when in doubt.*
    
The simplest usage of this package revolves around a single function: 
```julia
evolve!(p::AbstractParticle, bt::Vector{Obstacle}, total_time)
```
The function evolves a particle `p` inside a billiard table `bt` for a given amount of time `total_time`, 
while taking care of all the details internally. 

The first step is to define the billiard table `bt`, which is the system the particle `p` will propagate in. 
A billiard table is simply a collection (`Vector`) of `Obstacle`s. The most convenient way is to use 
one of the pre-defined billiard tables offered by the package. For example, let's create a periodic Sinai 
billiard with disk radius of 0.3 and with one side of length 2 and one of length 1:
```julia
bt = billiard_sinai(0.3, 2.0, 1.0; periodic=true)                                                    
```
(for more information about defining billiard tables see the official documentation here)

Afterwards, you want to create a particle inside that billiard system. 
For that, the function `randominside(bt::Vector{Obstacle})` is provided. 
This function returns a particle with random initial conditions inside the billiard table, 
while making sure that it is always in the allowed region of the billiard table.
```julia
p = randominside(bt)
```
Now you are ready to evolve this particle:
```julia
ct, poss, vels = evolve!(p, bt, 1000.0)
```
The return values of the `evolve!()` function need some brief explaining: As noted by the "!" at the end of the function, 
it changes its argument `p` (specifically, it updates almost all the fields of `p`).
Most importantly however, this function also returns the main output expected by a billiard
system. This output is a tuple of three vectors:
* `ct::Vector` : Collision times.
* `poss::Vector{SVector{2}}` : Positions during collisions.
* `vels:: Vector{SVector{2}})` : Velocities **exactly after** the collisions (i.e. reflections).

The time `t[i]` is the time necessary to reach state `poss[i], vels[i]` starting from the
state `poss[i-1], vels[i-1]`. That is why `t[1]` is always 0 since `poss[1], vels[1]` are
the initial conditions.

If this output is not convenient for you, the function `construct(t, poss, vels)` is provided, 
which constructs the (continuous) timeseries of the position and velocity, as well as the time-vector, when given the main output of `evolve!()`:
```julia
xt, yt, vxt, vyt, ts = construct(ct, poss, vels)
```
or, by taking advantage of the awesome ellipsis operator, you can do:
```julia
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0)...)
```
---
## Magnetic Propagation
The are only two differences between magnetic and straight propagation. 
Firstly, the particle type is not `Particle` anymore, but `MagneticParticle`. 
The latter has an extra field called `omega` which is the cyclic frequency of rotation 
(equivalently, the angular velocity). In order to create a `MagneticParticle` (without using the constructors), 
you simply provide this extra information to the `randominside()` function:
```julia
ω = 0.5
p = randominside(bt, ω)
typeof(p) # MagneticParticle
p.omega   # 0.5
```
To propagate the particle you use the same functions:
```julia
ct, poss, vels, ω = evolve!(p, bt, 1000.0)  #evolve for magnetic also returns ω
xt, yt, vxt, vyt, ts = construct(ct, poss, vels, ω, dt)
# or equivalently: 
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0)..., dt)
```
As you can see, the second difference is that the additional argument of the angular velocity must also be provided 
to the `construct()` function, in order for it to construct circular motion instead of straight motion between collisions. 
(Note: `evolve!()` returns 4 arguments for magnetic propagation, making the ellipsis syntax extremely useful!).

The final optional argument `dt` is the time-step at which the timeseries are constructed 
(since they are made up of sines and cosines).

---

## Ray-Splitting
No matter how complex ray-splitting processes you want, and irrespectively of
how many obstacles in the billiard table can perform ray-splitting, there is only
a single difference on the main function call:
The `evolve!()` function is supplemented with a fourth argument, called "ray_splitter":
```julia
ray_splitter::Dict{Int, Vector{Function}}
```
This argument is simply a dictionary which handles all ray-splitting processes in the billiard system.
It is a map of the Obstacle index within the billiard table to the
ray-splitting functions: (φ is the angle of incidence)
* T(φ, where, ω) : Transmission probability.
* θ(φ, where, ω) : Transmission (aka diffraction) angle.
* new_ω(ω, where) : Angular velocity after transmission.

Assuming you have defined a billiard table and a ray-splitter dictionary, the implementation is exactly the same as in the two previous cases: the ray-splitting dictionary is passed to `evolve!()` as a fourth argument.
```julia
ray_splitter = Dict(5 => [foo, bar, baz])
p = randominside(bt, 4.0)
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 100.0, ray_splitter)..., dt = 0.01)
```

For more information and instructions on defining the "ray_splitter" dictionary
please visit the "Ray-Splitting" tutorial here.

---
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
p = randominside(bt, ω)

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
