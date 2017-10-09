# Ray-Splitting
Ray-splitting is a semi-classical approach to the billiard system, giving a wave attribute to the ray traced by the particle.
Upon collision a particle may propagate through an obstacle (transmission & refraction) or be reflected. Following the mindset of this package, implementing a ray-splitting billiard requires only three very simple steps.

## Ray-Splitting Obstacles
First, an obstacle that supports ray-splitting is required to be present in your billiard table. The only new feature these obstacles have is an additional Boolean field called `pflag` (propagation flag). This field notes on which side of the obstacle the particle is currently propagating.

The normal vector as well as the distance from boundary change sign depending on the value of `pflag`. The obstacles `Antidot` and `SplitterWall` are the equivalents of disk and wall for ray-splitting.

Let's add an `Antidot` to a billiard table:

```julia
using DynamicalBilliards
bt = billiard_rectangle()
a = Antidot([0.5,0.5], 0.3)
push!(bt, a)
```

## Ray-Splitting Functions
Secondly, for each obstacle in your billiard table that you want to perform ray-splitting, you have to define 3 functions (φ = angle of incidence, ω = angular velocity *before* transmission):

1. T(φ, `pflag`, ω) : Transmission probability Τ depending on
   whether the particle is inside or outside the obstacle (`pflag`) and optionally depending on ω.
2. θ(φ, `pflag`, ω) : Transmission (aka refraction) angle θ
   depending on whether the particle is inside or outside the obstacle (`pflag`) and optionally depending on ω.
3. ω_new(ω, `pflag`) : Angular velocity after transmission.

The above three functions use the **same convention**: the argument `pflag` is the one the Obstacle has **before transmission**. For example, if a particle is outside a disk (let `pflag = true` here) and is transmitted inside the disk (`pflag` becomes `false` here), then all three functions will be given their second argument (the Boolean one) as `true`!

## Ray-Splitter Dictionary
To pass the information of the aforementioned functions into the main API a dictionary is required:
```julia
raysplitter::Dict{Int, Any}
```
This dictionary is a map of the **obstacle index** within the billiard table to a **container of the ray-splitting functions**. This container could be a `Vector` or a `Tuple` and the later is the suggested version.

For example, if we wanted to allocate ray-splitting functions for the **5th** obstacle in our billiard table, which could be e.g. an `Antidot`, we would write something like:
```julia
sa = (θ, pflag, ω) -> pflag ? 2θ : 0.5θ  # refraction (scatter) angle
T = (θ, pflag, ω) -> begin   # Transmission probability
  if pflag
    abs(θ) < π/4 ? 0.5exp(-(θ)^2/2(π/8)^2) : 0.0
  else
    0.75*exp(-(θ)^2/2(π/4)^2)
  end
end
newo = (ω, bool) -> bool ? -0.5ω : -2ω   # new angular velocity
raysplitter = Dict(5 => (T, sa, newo))  # Index maps to container of Functions
```

!!! note "Order of Arguments"
    The functions **must accept the specific number of arguments shown in the previous section** even if some are not used. Also, the functions must be given **in the specific order: [1. transmission probability, 2. refraction angle, 3. new ω]** in the vector passed to the dictionary.

The next step is very simple: the `raysplitter` dictionary is directly passed into `evolve!()` as a fourth argument.
Using the billiard table we defined previously, where its 5th element is a ray-splitting `Antidot`, we now do:
```julia
bt = billiard_rectangle()
a = Antidot([0.5, 0.5], 0.25)
push!(bt, a)
ω = 1.25
p = randominside(bt, ω)
dt = 0.05
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 25.0, raysplitter)..., dt)
using PyPlot
DynamicalBilliards.enableplotting()
plot_billiard(bt)
plot(xt, yt)
```
which should give a result similar to this:

![Ray-splitting billiard](http://i.imgur.com/UfGQfOm.png)

A final difference to be noted is that in the case of ray-splitting with magnetic fields, the fourth value returned by `evolve!()` is not a number, but a vector of angular velocities `omegas`. The value `omegas[i]` is the angular velocity the particle has while propagating from state `pos[i], vel[i]` to state `pos[i+1], vel[i+1]`. The `construct()` function takes this into account.

### No field `pflag` error

If you ever encounter the error `ERROR: type SomeObstacleType has no field pflag` this means that the index provided by your ray-splitting dictionary points to an object that does not support ray-splitting. Use the functions:
```julia
acceptable_raysplitter(raysplitter, bt)
supports_raysplitting(obst::Obstacle)
```
to find out what you did wrong. Most likely, the index you supplied was incorrect, i.e. the index could be `5` instead of `4`.

## Example Animation
In the [examples page], you can find the code for the following animation, which
includes ray-splitting:

![Ray-splitter animation](http://i.imgur.com/89s0fon.gif)

## Physics
The condition for transmission is simply: `T(φ, pflag, ω) > rand()`. If it returns `true`, transmission (i.e. ray-splitting) will happen. Otherwise just specular reflection will take place. A more detailed discussion is on the ray-splitting section of the
[Physics](/physics#ray-splitting-functions) page.
