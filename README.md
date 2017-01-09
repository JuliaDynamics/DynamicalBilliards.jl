# DynamicalBilliards.jl
A Julia package for dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction.

The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:
* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard. 
* **Ray-splitting Billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiard systems that have already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : functions for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more. Also includes animation related content.

NOTICE: This package does not support collision between particles. All particles are considered point particles for all simulations offered by `DynamicalBilliards.jl`.

Documentation: link.

*(all exported names of DynamicalBilliards.jl have detailed docstrings. Use `?` when in doubt)*

## Basic Usage
The simplest usage of this package revolves around a single function: 
```
evolve!(p::AbstractParticle, bt::Vector{Obstacle}, total_time)
```
The function evolves a particle `p` inside a billiard table `bt` for a given amount of time `total_time`, while taking care of all the details internally. 

The first step is to define the billiard table `bt`, which is the system the particle `p` will propagate in. A billiard table is simply a collection (`Vector`) of `Obstacle`s. The most convenient way is to use one of the pre-defined billiard tables offered by the package. For example, let's create a periodic Sinai billiard with disk radius of 0.3 and with one side of length 2 and one of length 1:
```
bt = billiard_sinai_periodic(0.3, 2.0, 1.0)                                                    
```
(for more information about defining billiard tables see the official documentation here)

Afterwards, you would want to create a particle inside that billiard system. For that, the function `randominside(bt::Vector{Obstacle})` is provided. This function returns a particle with random initial conditions inside the billiard table, while making sure that it is always in the allowed region of the billiard table.
```
p = randominside(bt)
```
Now you are ready to evolve this particle:
```
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
```
xt, yt, vxt, vyt, ts = construct(t, poss, vels)
```
or, if you want to use the fancy ellipsis operator, you can do:
```
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 1000.0)...)
```

## Visualizing

## Magnetic Propagation
