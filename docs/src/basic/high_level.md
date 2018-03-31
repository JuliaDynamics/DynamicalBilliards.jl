# High Level API

`DynamicalBilliards` was created with easy-of-use as its main cornerstone.
With 3 simple steps, the user can get the output of the propagation of a particle in a billiard.

In general, the workflow of `DynamicalBilliards` follows these simple steps:
1. Create a billiard.
2. Create a particle inside that billiard.
3. Get the output you want by using one of the high level functions.

Adding more complexity in your billiard does not add complexity in your code. For example, to implement a ray-splitting billiard
you only need to define one additional variable, a dictionary `Dict{Int, Vector{Function}}`.

After reading through this basic usage page,
you will be able to use all aspects of `DynamicalBilliards.jl` with minimal effort.

!!! tip "Visualizations"
    Visualizing the billiards, particles and their motion is one of the most important parts of the `DynamicalBilliards`. It is not discussed in this page however, but rather in the [Visualizing](/visualizing) page.

---
## Billiard
A [`Billiard`](@ref) is simply a collection of [`Obstacle`](@ref) subtypes. Particles are propagating inside a `Billiard`, bouncing from obstacle to obstacle while having constant velocity in-between.

There is a [tutorial](/tutorials/billiard_table) on how to create your own billiard. In addition, there are many pre-defined billiards that can be found in the [Standard Billiards Library](#standard-billiards-library) section. That is why knowing how to construct a [`Billiard`](@ref) is not important at this point.

In this page we will be using the Bunimovich billiard as an example:
```julia
using DynamicalBilliards
bt = billiard_bunimovich() # using default arguments
```
```
Billiard{Float64} with 4 obstacles:
  Bottom wall
  Right semicircle
  Top wall
  Left semicircle
```

## Particles
A "particle" is that thingy that moves around in the billiard. It always moves with velocity measure 1, by convention.

Currently there are two types of particles:

* `Particle`, which propagates as a straight line.
* `MagneticParticle`, which propagates as a circle instead of a line (similar to electrons in a perpendicular magnetic field).

There are two ways to create a particle. The first one is to provide the
constructor with some initial conditions:
```julia
x0 = rand(); y0 = rand();
φ0 = 2π*rand() # an angle is enough
p = Particle(x0, y0, φ0)
```
```
Particle{Float64}
position: [0.324647, 0.142048]
velocity: [0.573107, 0.81948]
```
To create a `MagneticParticle` simply provide the constructor with one more number,
the angular velocity:
```julia
ω = 0.5
mp = MagneticParticle(x0, y0, φ0, ω)
```
```
MagneticParticle{Float64}
position: [0.324647, 0.142048]
velocity: [0.573107, 0.81948]
ang. velocity: 0.5
```


!!! faq "Why the `{Float64}` ?"
    When creating a billiard or a particle, the object is printed with `{Float64}` at the end. This shows what type of numbers are used for *all* numerical operations. If you are curious you can learn more about it in the [numerical precision page](/physics/#numerical-precision).

!!! danger "Particles must be inside the Billiard!"
    Keep in mind that the particle must be initialized **inside a billiard** for any functionality to work properly and make sense. If you are not sure what we mean by that, then you should check out the [low-level API page](LINKME).

## Random initial conditions

If you have a `Billiard` which is not a rectangle, creating many random initial conditions inside it can be a pain. Fortunately, the second way to create a particle is to use the following function:
```@docs
randominside
```
---

For example:
```julia
p = randominside(bt)
```
```
Particle{Float64}
position: [0.274096, 0.612643]
velocity: [-0.178995, 0.98385]
```
and
```julia
mp = randominside(bt, ω)
```
```
MagneticParticle{Float64}
position: [0.267054, 0.631786]
velocity: [-0.500439, -0.865772]
ang. velocity: 0.5
```
`randominside` always creates particles with same number type as the billiard.

## `evolve` & `construct`
Now that we have created a billiard and a particle inside, we want to evolve it!
There is a simple function for that, called `evolve!` (or `evolve` if you don't want to mutate the particle):
```@docs
evolve!
```
---
Forget the ray-splitting part for now, we'll talk a bit more about it later on (see [Ray-Splitting](#ray-splitting)).

Let's see an example:
```julia
ct, poss, vels = evolve(p, bt, 100)
for i in 1:5
  println(round(ct[i], 3), "  ", poss[i], "  ", vels[i])
end
```
```
0.0    [0.274096, 0.612643]  [-0.178995, 0.98385]
0.394  [0.203623, 1.0]  [-0.178995, -0.98385]
1.016  [0.0216906, 0.0]  [-0.178995, 0.98385]
0.991  [-0.155718, 0.975134]  [0.438064, -0.898944]
1.085  [0.319474, 0.0]  [0.438064, 0.898944]
```

Similarly, for magnetic propagation
```julia
ct, poss, vels, ω = evolve(mp, bt, 100)
for i in 1:10
  println(round(ct[i], 3), "  ", poss[i], "  ", vels[i])
end
```
```
0.0    [0.267054, 0.631786]  [-0.500439, -0.865772]
0.677  [0.0329496, 2.94431e-13]  [-0.184546, 0.982824]
0.947  [-0.351508, 0.855587]  [0.783478, -0.621419]
1.959  [1.4997, 0.517318]  [-0.971386, 0.237507]
1.905  [-0.283005, 0.0878011]  [0.338364, 0.941015]
```

The above return types are helpful in some applications. In other applications however
one prefers to have the time series of the individual variables. For this,
the `construct` function is used:
```@docs
construct
```
---
The function is especially useful when one wants immediately the timeseries instead
of the output of `evolve!`. Because of the simple syntax
```julia
xt, yt, vxt, vyt, t = construct(evolve(p, bt, 100)...)

# print as a matrix:
hcat(xt, yt, vxt, vyt, t)[1:5, :]
```
```
5×5 Array{Float64,2}:
  0.274096   0.612643  -0.178995   0.98385   0.0     
  0.203623   1.0       -0.178995  -0.98385   0.393715
  0.0216906  0.0       -0.178995   0.98385   1.41013
 -0.155718   0.975134   0.438064  -0.898944  2.40127
  0.319474   0.0        0.438064   0.898944  3.48603
```

This nicely reveals why in the case of magnetic propagation `evolve!` also returns
the angular velocity. So that it is possible to do the same process for magnetic
propagation as well (plus, it is also useful in ray-splitting).
```julia
# evolve the magnetic particle instead:
xt, yt, vxt, vyt, t = construct(evolve(mp, bt, 100)...)

# print as a matrix:
hcat(xt, yt, vxt, vyt, t)[1:5, :]
```
```
5×5 Array{Float64,2}:
 0.267054  0.631786  -0.500439  -0.865772  0.0
 0.262071  0.623116  -0.496104  -0.868263  0.01
 0.257132  0.614421  -0.491756  -0.870733  0.02
 0.252236  0.605701  -0.487397  -0.873181  0.03
 0.247384  0.596957  -0.483025  -0.875607  0.04
```


!!! note "Type of `t`"
    Remember that the behavior of `evolve!` depends on the type of the third argument,
    which represents "total amount". If it is `AbstractFloat`, it represents total amount of time, but if it is `Int` it represents total number of collisions.


## Poincaré Sections
Poincaré sections (also known as boundary maps) can be obtained with the high
level function
```@docs
poincaresection
```
---
For example, take a look at boundary maps of the mushroom billiard, which is known to have a mixed phase space:
```julia
using DynamicalBilliards

bt = billiard_mushroom()

n = 100 # how many particles to create

ξς, φς, ις = poincaresection(bt, 10000, n)

using PyPlot # enables plot_psos function

colors = ["C$(rand(1:9))" for i in 1:n] # random colors

figure()
plot_psos(ξς, φς, ις, color = colors)
```
![PSOS1](https://i.imgur.com/RO9UZa9.png)

And of course similarly for magnetic fields
```julia
ξς, φς, ις = poincaresection(bt, 10000, n, 1.0) # angular velocity last argument
figure()
plot_psos(ξς, φς, ις, color = colors)
```
![PSOS2](https://i.imgur.com/YoW1FVD.png)

## Escape Times
It is very easy to create your own function that calculates an "escape time": the time until the particle leaves the billiard by meeting a specified condition. There is also a high-level function for this though:
```@docs
escapetime
```
!!! tip "Creating a door"
    To create a "door" simply visit the [library page](library) to learn more about the individual obstacle types (specifically [`FiniteWall`](@ref)). To be able to
    combine them into a [`Billiard`](@ref) you should also check out the tutorial on [defining your own billiard](tutorials/billiard_table).

For example, the default implementation of the mushroom billiard has a "door" at the
bottom of the stem. Thus,
```julia
bt = billiard_mushroom()
et = zeros(100)
for i ∈ 1:100
    p = randominside(bt)
    et[i] = escapetime(p, bt, 10000)
end
println("Out of 100 particles, $(count(x-> x != Inf, et)) escaped")
println("Mean escape time was $(mean(et[et .!= Inf]))")
```
```
Out of 100 particles, 11 escaped
Mean escape time was 4.436943929428599
```
Of course, `escapetime` works with `MagneticParticle` as well
```julia
escapetime(randominside(bt, 1.0), bt, 10000)
```
```
87.17643414941281
```



## Ray-Splitting
During ray-splitting a particle may propagate
through an obstacle given arbitrary transmission and refraction
laws. This is also known as a "semiclassical billiard".

No matter how complex ray-splitting processes you want, and irrespectively of
how many obstacles in the billiard table can perform ray-splitting, there is only
a single difference on the main function call:
The `evolve!()` function is supplemented with a fourth argument:
```julia
ray_splitter::Dict{Int, Any}
```
This argument is simply a dictionary which handles all ray-splitting processes in the billiard system.

Given an obstacle index within the billiard table it returns a container of the
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

For more information and instructions on defining the `ray_splitter` dictionary visit the [Ray-Splitting tutorial](/tutorials/ray-splitting).


## Standard Billiards Library
```@docs
billiard_rectangle
billiard_sinai
billiard_bunimovich
billiard_mushroom
billiard_lorentz
billiard_polygon
billiard_hexagonal_sinai
billiard_raysplitting_showcase
```
