#Basic Usage

`DynamicalBilliards.jl` was created with easy-of-use as its main cornerstone.
With 3 simple steps, the user can fully initalize, evolve, and get the output of the propagation of a particle in a billiard system.

In general, the workflow of `DynamicalBilliards.jl` follows these simple steps:
1. Create a billiard table, a `Vector{Obstacle}`.
2. Create a particle inside that billiard table.
3. Get the output by evolving the particle.

Adding more complexity in your billiard table does not add complexity in your code. For example, to implement a ray-splitting billiard
you only need to define one additional variable, a dictionary `Dict{Int, Vector{Function}}`. After reading through this basic usage page,
you will be able to use all aspects of `DynamicalBilliards.jl` with minimal effort.

---

## Straight Propagation

The usage of this package revolves around a single function:
```julia
evolve!(p::AbstractParticle, bt::Vector{Obstacle}, t::Union{Int, Float})
```
which evolves a particle `p` inside a billiard table `bt`. If the given `t` is of type `AbstractFloat`, the evolution happens for `t` amount of time. If however `t` is of type `Int`, the evolution happens for `t` number of collisions (other types are not supported).

The first step is to define the billiard table `bt`, which is the system the particle `p` will propagate in.
A billiard table is simply a collection (`Vector`) of `Obstacle`s. The most convenient way is to use
one of the pre-defined billiard tables offered by the package. For example, let's create a periodic Sinai
billiard with disk radius of 0.3 and with one side of length 2 and one of length 1:
```julia
using DynamicalBilliards
bt = billiard_sinai(0.3, 2.0, 1.0; setting = "periodic")
```
For more information about defining billiard tables visit the [tutorial on defining your own billiard table](/tutorials/billiard_table)). You should definitely look up that page
if you want to customly define a table instead of using the [predefined ones](basic/library/#standard-billiards).

Afterwards, you want to create a particle inside that billiard system.
For that, the function `randominside(bt::Vector{Obstacle})` is provided, which returns a particle with random initial conditions inside the billiard table.
```julia
p = randominside(bt)
```
If you want to specify the initial conditions yourself, simply pass them to the `Particle` constructor, like `p = Particle(x0, y0, φ0)`.
Now you are ready to evolve this particle:
```julia
t = 1000.0 # subtype of AbstractFloat
ct, poss, vels = evolve!(p, bt, 1000.0)
```

!!! note "Type of `t`"
    The behavior of `evolve!` depends on the type of the third argument,
    which represents total amount. If it is `AbstractFloat`, it represents total amount of time, but if it is `Int` it represents total numnber of collisions.

The return values of the `evolve!()` function need some brief explaining: As noted by the "!" at the end of the function, it changes its argument `p`.
Most importantly however, this function also returns the main output expected by a billiard system. This output is a tuple of three vectors:
* `ct::Vector{T}` : Collision times.
* `poss::Vector{SVector{2,T}}` : Positions during collisions.
* `vels::Vector{SVector{2,T}}` : Velocities **exactly after** the collisions (e.g. after reflection).

with `T<:AbstractFloat` (more on that on the [Numerical Precision](/physics/#numerical-precision) section). The time `ct[i]` is the time necessary to reach state `poss[i], vels[i]` starting from the
state `poss[i-1], vels[i-1]`. That is why `ct[1]` is always 0 since `poss[1], vels[1]` are
the initial conditions.

If this output is not convenient for you, the function `construct(ct, poss, vels)` is provided,
which constructs the (continuous) timeseries of the position and velocity, as well as the time-vector, when given the main output of `evolve!()`:
```julia
xt, yt, vxt, vyt, ts = construct(ct, poss, vels)
```
or, by taking advantage of the awesome ellipsis operator, you can do:
```julia
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 100.0)...)
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
typeof(p) # MagneticParticle{Float64}
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

The final *optional argument* `dt` is the time-step at which the timeseries are constructed
(since they are made up of sines and cosines).

---

## Ray-Splitting

No matter how complex ray-splitting processes you want, and irrespectively of
how many obstacles in the billiard table can perform ray-splitting, there is only
a single difference on the main function call:
The `evolve!()` function is supplemented with a fourth argument, called "ray_splitter":
```julia
ray_splitter::Dict{Int, Any}
```
This argument is simply a dictionary which handles all ray-splitting processes in the billiard system.
It is a map of the Obstacle index within the billiard table a container of the
ray-splitting functions: (φ is the angle of incidence)
* T(φ, pflag, ω) : Transmission probability.
* θ(φ, pflag, ω) : Transmission (aka diffraction) angle.
* new_ω(ω, pflag) : Angular velocity after transmission.

Assuming you have defined a billiard table and a ray-splitter dictionary, the implementation is exactly the same as in the two previous cases: the ray-splitting dictionary is passed to `evolve!()` as a fourth argument.
```julia
ray_splitter = Dict(5 => (foo, bar, baz))
p = randominside(bt, 4.0)
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 100.0, ray_splitter)..., 0.01)
```

For more information and instructions on defining the "ray_splitter" dictionary visit the [Ray-Splitting tutorial here](/tutorials/ray-splitting).

---

## Visualizing

The functions `plot_obstacle(obst::Obstacle; kwargs...)`, `plot_billiard(bt::Vector{Obstacle})` and `plot_particle(p::AbstractParticle; kwargs...)` are provided in order to plot the respective elements **on the current PyPlot figure**. The `kwargs...` are keywords passed directly into `PyPlot`'s constructors (like e.g. `linewidth = 2.0`).

[The tutorial on visualizing](/tutorials/visualizing) has step-by-step descriptions on how to handle all plotting offered by `DynamicalBilliards.jl`.

### Introduction animation

The example .gif shown in the introduction, was generated simply with the code:
```julia
using DynamicalBilliards, PyPlot

bt = billiard_rectangle(1.5, 1.0)
d1 = Disk([0.45, 0.6], 0.3, "Upper-left Disk")
d2 = Disk([1.1, 0.3], 0.15, "Lower-right Disk")
d3 = Disk([1.2, 0.8], 0.1, "Small Disk")
w1 = FiniteWall([0.0, 0.4], [0.6,0.0], [0.4,0.6], "Diagonal")
push!(bt, d1, d2, d3, w1)
ω = 2.0
p = randominside(bt, ω)

DynamicalBilliards.enableplotting()
plot_billiard(bt)
axis("off")
tight_layout()

sname = "C:\\***\\example"
animate_evolution(p, bt, 200;
col_to_plot = 4, savefigs = true, savename = sname)
```
Afterwards the outputed .png files where merged into a single .gif externally using for example the website gifmaker.me.
