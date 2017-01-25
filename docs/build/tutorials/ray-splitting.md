
<a id='Ray-Splitting-1'></a>

# Ray-Splitting


Ray-splitting is a semiclassical approach to the billiard system, giving a wave attribute to the ray traced by the particle. Upon collision a particle may propagate through an obstacle (transmission & refraction) or be reflected. Following the mindset of this package, implementing a ray-splitting billiard requires only three very simple steps.


<a id='Ray-Splitting-Obstacles-1'></a>

## Ray-Splitting Obstacles


First, an obstacle that supports ray-splitting is required to be present in your billiard table. The only new feature these obstacles have is an additional Boolean field called `where`. This field notes on which side of the obstacle the particle is currently propagating *(if you are wondering how a distance can change sign, see the tutorial on Defining your own Obstacles)*. The normal vector as well as the distance from boundary change sign depending on the value of `where`. The obstacles `Antidot` and `SplitterWall` are the equivalents of disk and wall for ray-splitting. To make your own defined obstacle support ray-splitting, visit this tutorial.


> There is a simple reason for having extra Types to support ray-splitting: non ray-splitting Types  always perform 2 less operations in their innermost loops, saving a bit of time. Also, ray-splitting obstacles are defined as `type` instead of as `immutable`.



Let's add an `Antidot` to a billiard table:


```julia
using DynamicalBilliards
bt = billiard_rectangle()
a = Antidot([0.5,0.5], 0.3)
push!(bt, a)
```


<a id='Ray-Splitting-Functions-1'></a>

## Ray-Splitting Functions


Secondly, for each obstacle in your billiard table that will perform ray-splitting, you have to define 3 functions. Notice that not every obstacle that supports ray-splitting actually has to perform it; it is up to the user. Those 3 functions are the following:


1. T(φ, `where`, ω) : Takes as input the angle of incidence φ and returns the transmission probability Τ depending on  whether the particle is inside or outside the obstacle (`where`) and optionally depending on ω.  This function should be an even function with respect to φ.
2. θ(φ, `where`, ω) : Takes as input the angle of incidence 	φ and returns the the transmission (aka refraction)  angle θ  depending on whether the particle is inside or outside the obstacle (`where`) and optionally depending on ω.  This function should be an odd function with respect to φ.
3. ω_new(ω, `where`) : Angular velocity after transmission.


The above three functions use the **same convention**: the argument `where` is the one the Obstacle has **before transmission**. For example, if a particle is outside a disk (let `where = true` here) and is transmitted inside the disk (`where` becomes `false` here), then all three functions will be given their second argument (the boolean one) as `true`!


<a id='Ray-Splitter-Dictionary-1'></a>

## Ray-Splitter Dictionary


To pass the information of the aforementioned functions into the main API (`evolve!()`) a dictionary is required, which we will call "raysplitter": `raysplitter::Dict{Int, Vector{Function}}`. The keys are integers and the values are vectors of functions.  This dictionary is a map of the obstacle index within the billiard table to the ray-splitting functions. For example, if we wanted to allocate ray-splitting functions for the 5th obstacle in our billiard table, which could be e.g. an `Antidot`, we would write something like:


```julia
sa = (θ, where, ω) -> where ? 2θ : 0.5θ
T = (θ, where, ω) -> begin
  if where
    abs(θ) < π/4 ? 0.5exp(-(θ)^2/2(π/8)^2) : 0.0
  else
    0.75*exp(-(θ)^2/2(π/4)^2)
  end
end
newo = (ω, bool) -> bool ? -0.5ω : -2ω
raysplitter = Dict(5 => [T, sa, newo])
```


Notice the following two very important points: The functions **must accept the specific number of arguments shown in the previous section** even if some are not used. Also, the functions must be given **in the specific order: 1. transmission probability, 2. refraction angle, 3. new ω** in the vector passed to the dictionary.


The next step is very simple: the `raysplitter` dictionary is directly passed into `evolve!()` as a fourth argument.  Using the billiard table we defined previously, where its 5th element is a ray-splitting `Antidot`, we now do:


```julia
ω = 4.0
p = randominside(bt, ω)
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 100.0, raysplitter)..., dt = 0.05)
plot_billiard(bt)
plot(xt, yt)
```


and everything works like a charm! A final difference to be noted is that in the case of ray-splitting with magnetic fields, the fourth value returned by `evolve!()` is not a number, but a vector of angular velocities `omegas`. The value `omegas[i]` is the angular velocity the particle has while propagating from state `pos[i], vel[i]` to state `pos[i+1], vel[i+1]`. The `construct()` function takes this into account.


<a id='Physics-1'></a>

## Physics


The condition for transmission is simply: `T(φ, where, ω) > rand()`. If it returns `true`, transmission (i.e. ray-splitting) will happen. Otherwise just specular reflection will take place.


The three key functions given to the `raysplitter` dictionary must have some properties in order to have physical meaning, like for example that the scattering probability function is even towards φ. One of these properties is absolutely **mandatory** for this package to work properly. This is the property of total internal reflection, i.e. if the refraction angle is calculated to be greater/equal than π/2, no transmission can happen. **This condition is not assured internally** and thefore you must be sure that your transmission probability function satisfies it. In the above example, the function `T` makes sure to return 0 in that case.


In order to test if the `raysplitter` dictionary you have defined has physical meaning, the function `isphysical()` is provided. Its documentation string has all the details one should know:


> ```julia
> isphysical(raysplitter::Dict{Int, Vector{Function}}; only_mandatory = false)
> ```
>
> Return `true` if the given ray-splitting dictionary has properties of the physical world.
>
> Specifically, check if (φ is the incidence angle):
>
>   * Critical angle means total reflection: If θ(φ) ≥ π/2 then T(φ) = 0
>   * Transmission probability is even function: T(φ) ≈ T(-φ)
>   * Refraction angle is odd function: θ(φ) ≈ -θ(-φ)
>   * Ray reversal is true: θ(θ(φ, where, ω), !where, ω) ≈ φ
>   * Magnetic conservation is true: (ω_new(ω_new(ω, where), !where) ≈ ω
>
> The first property is mandatory and must hold for correct propagation. They keyword `only_mandatory` notes whether the rest of the properties should be tested or not.



<a id='Velocity-measure-1'></a>

### Velocity measure


The `MagneticParticle` is assumed to always have a velocity vector of measure 1. This simplyfies the formulas used internally to a significant amount. During ray-splitting, the particle may be in areas with different fields, resulting in different angular velocities (result of the "ω_new" function).


Physically, in such a situation, the velocity measure of the particle will also change. This is not implemented internally by `DynamicalBilliards`  during the evolution of the particle. However it is very easy to change this after the particle has evolved.

