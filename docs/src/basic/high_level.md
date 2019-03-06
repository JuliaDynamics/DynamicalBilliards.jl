# High Level API

`DynamicalBilliards` was created with ease-of-use as its main cornerstone.
With 3 simple steps, the user can get the output of the propagation of a particle in a billiard.

In general, the workflow of `DynamicalBilliards` follows these simple steps:
1. Create a billiard.
2. Create particles inside that billiard.
3. Get the output you want by using one of the high level functions.

Adding more complexity in your billiard does not add complexity in your code.
For example, to implement a ray-splitting billiard you only need to define one additional variable, a [`RaySplitter`](@ref) and pass it to the high level functions.

After reading through this page, you will be able to use almost all aspects of `DynamicalBilliards` with minimal effort.

!!! tip "Visualizations"
    Visualizing the billiards, particles and their motion is one of the most important parts of the `DynamicalBilliards`. It is not discussed in this page however, but rather in the [Visualizing](/visualizing) page.

---
## Billiard
A [`Billiard`](@ref) is simply a collection of [`Obstacle`](@ref) subtypes. Particles are propagating inside a `Billiard`, bouncing from obstacle to obstacle while having constant velocity in-between.

There is a [tutorial](/tutorials/billiard_table) on how to create your own billiard. In addition, there are many pre-defined billiards that can be found in the [Standard Billiards Library](#standard-billiards-library) section. That is why knowing how to construct a [`Billiard`](@ref) is not important at this point.

In this page we will be using the Bunimovich billiard as an example:
```@example 2
using DynamicalBilliards
bd = billiard_bunimovich()
```

## Particles
A "particle" is that thingy that moves around in the billiard. It always moves with velocity of measure 1, by convention.

Currently there are two types of particles:

* [`Particle`](@ref), which propagates as a straight line.
* [`MagneticParticle`](@ref), which propagates as a circle instead of a line (similar to electrons in a perpendicular magnetic field).

There are two ways to create a particle. The first one is to provide the
constructor with some initial conditions:
```@example 2
x0 = rand(); y0 = rand();
φ0 = 2π*rand()
p = Particle(x0, y0, φ0)
```

To create a `MagneticParticle` simply provide the constructor with one more number,
the angular velocity:
```@example 2
ω = 0.5
mp = MagneticParticle(x0, y0, φ0, ω)
```


!!! faq "Why the `{Float64}` ?"
    When creating a billiard or a particle, the object is printed with `{Float64}` at the end. This shows what type of numbers are used for *all* numerical operations. If you are curious you can learn more about it in the [numerical precision page](/low_level/#numerical-precision).

!!! danger "Particles must be inside the Billiard!"
    Keep in mind that the particle must be initialized **inside a billiard** for any functionality to work properly and make sense. If you are not sure what we mean by that, then you should check out the [Internals page](low_level).

## Random initial conditions

If you have a `Billiard` which is not a rectangle, creating many random initial conditions inside it can be a pain. Fortunately, the second way to create a particle is to use the following function:
```@docs
randominside
```
---

For example:
```@example 2
p = randominside(bd)
```

and
```@example 2
mp = randominside(bd, ω)
```

`randominside` always creates particles with same number type as the billiard.

## `evolve` & `timeseries`
Now that we have created a billiard and a particle inside, we want to evolve it!
There is a simple function for that, called `evolve!` (or `evolve` if you don't want to mutate the particle), which returns the time, position and velocities at the collision points:
```@docs
evolve!
```
---
Forget the ray-splitting part for now (see [Ray-Splitting](/ray-splitting)).

Let's see an example:
```@example 2
ct, poss, vels = evolve(p, bd, 100)
for i in 1:5
  println(round(ct[i], digits=3), "  ", poss[i], "  ", vels[i])
end
```

Similarly, for magnetic propagation
```@example 2
ct, poss, vels, ω = evolve(mp, bd, 100)
for i in 1:10
  println(round(ct[i], digits=3), "  ", poss[i], "  ", vels[i])
end
```

The above return types are helpful in some applications.
In other applications however one prefers to have the time series of the individual variables.
For this, the `timeseries` function is used:
```@docs
timeseries!
```
---
For example:
```@example 2
xt, yt, vxt, vyt, t = timeseries(p, bd, 100)

# print as a matrix:
hcat(xt, yt, vxt, vyt, t)[1:5, :]
```

Same story for magnetic particles:
```@example 2
# evolve the magnetic particle instead:
xt, yt, vxt, vyt, t = timeseries(mp, bd, 100)

# print as a matrix:
hcat(xt, yt, vxt, vyt, t)[1:5, :]
```

Sometimes we may need information about which obstacles a particle visited, in which sequence, and when. For this we have the following function:
```@docs
visited_obstacles!
```

!!! note "Type of `t`"
    Remember that the behavior of time evolution depends on the type of the `t` argument, which represents "total amount". If it is `AbstractFloat`, it represents total amount of time, but if it is `Int` it represents total number of collisions.



## Poincaré Sections
```@docs
psos
```
---
For example, the surface of section in the periodic Sinai billiard with magnetic field
reveals the mixed nature of the phase-space:
```@example psos
using DynamicalBilliards, PyPlot
t = 100; r = 0.15
bd = billiard_sinai(r, setting = "periodic")

# the direction of the normal vector is IMPORTANT!!!
# (always keep in mind that ω > 0  means counter-clockwise rotation!)
plane = InfiniteWall([0.5, 0.0], [0.5, 1.0], [-1.0, 0.0])

posvector, velvector = psos(bd, plane, t, 1000, 2.0)
c(a) = length(a) == 1 ? "C1" : "C0"

figure()
for i in 1:length(posvector)
    poss = posvector[i] # vector of positions
    vels = velvector[i] # vector of velocities at the section
    L = length(poss)
    if L > 0
        #plot y vs vy
        y = [a[2] for a in poss]
        vy = [a[2] for a in vels]

        plot(y, vy, ls = "None", color = c(y), ms = 2.0, alpha = 0.75, marker = "o")
    end
end
xlabel("\$y\$"); ylabel("\$v_y\$")
savefig("psos.png"); nothing # hide
```
![](psos.png)

!!! note "`psos` operates on the unit cell"
    The `psos` function always calculates the crossings *within* the unit cell of
    a periodic billiard. This means that no information about the "actual" position
    of the particle is stored, everything is modulo the unit cell.

    This can be seen very well in the above example, where there aren't any entries
    in the region `0.5 - r ≤ y ≤ 0.5 + r`.

Of course it is very easy to "re-normalize" the result such that it is coherent.
The only change we have to do is simply replace the line `y = [a[2] for a in poss]`
with
```julia
y = [a[2] < 0.5 ? a[2] + 1 : a[2]  for a in poss]
```

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
```@example 2
using Statistics
bd = billiard_mushroom()
et = zeros(100)
for i ∈ 1:100
    particle = randominside(bd)
    et[i] = escapetime(particle, bd, 10000)
end
println("Out of 100 particles, $(count(x-> x != Inf, et)) escaped")
println("Mean escape time was $(mean(et[et .!= Inf]))")
```

Of course, `escapetime` works with `MagneticParticle` as well
```@example 2
escapetime(randominside(bd, 1.0), bd, 10000)
```


## Mean Collision Times
Because the computation of a mean collision time (average time between collisions
in a billiard) is often a useful quantity, the following function computes it
```@docs
meancollisiontime
```
For example,
```@example 2
bd = billiard_sinai()
meancollisiontime(randominside(bd), bd, 10000.0)
```

## Parallelization
```@docs
parallelize
```
---
Here are some examples
```@example 2
bd = billiard_stadium()
particles = [randominside(bd) for i in 1:1000]
parallelize(meancollisiontime, bd, 1000, particles)
```

```@example 2
parallelize(lyapunovspectrum, bd, 1000, particles)
```


## It's all about bounce!
The main propagation algorithm used by `DynamicalBilliards` is bundled in the following well-behaving function:
```@docs
bounce!
```
---
`bounce!` is the function used internally by all high-level functions, like [`evolve!`](@ref), [`boundarymap`](@ref), [`escapetime`](@ref), etc.

This is the function a user should use if they want to calculate other things besides what is already available in the high level API.


## Standard Billiards Library
!!! tip "You can also use keywords!"
    All standard billiards have a function version that accepts keyword arguments instead of positional arguments, for ease of use.

```@docs
billiard_rectangle
billiard_sinai
billiard_bunimovich
billiard_mushroom
billiard_polygon
billiard_hexagonal_sinai
billiard_raysplitting_showcase
billiard_logo
billiard_iris
```

## Particle types
```@docs
Particle
MagneticParticle
```
