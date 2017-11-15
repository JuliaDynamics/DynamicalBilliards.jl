---
title: 'DynamicalBilliards.jl: An easy-to-use, modular and extendable Julia package for Dynamical Billiard systems in two dimensions.'
tags:
  - billiards
  - physics
  - chaos
  - escape
  - obstacle
authors:
 - name: George Datseris
   orcid: 0000-0002-6427-2385
   affiliation: 1
affiliations:
 - name: Max Planck Institute for Dynamics and Self-Organization
   index: 1
date: 26 October 2017
bibliography: none
---

![Logo of DynamicalBilliards.jl](DynBil_logo---Animated.gif)

# Summary

[DynamicalBilliards.jl](http://orcid.org/0000-0002-6427-2385) is a package about
two-dimensional (dynamical) billiard systems written in its entirety in Julia. It is easy to use, and easy to
be extended. It is accompanied by a detailed documentation archive with a lot of tutorials, example code as well as information about the physical algorithms at play.

Documentation:

* [stable](https://juliadynamics.github.io/DynamicalBilliards.jl/stable/): https://juliadynamics.github.io/DynamicalBilliards.jl/stable/
* [latest](https://juliadynamics.github.io/DynamicalBilliards.jl/latest/): https://juliadynamics.github.io/DynamicalBilliards.jl/latest/

Repository page:

* [GitHub](https://github.com/JuliaDynamics/DynamicalBilliards.jl): https://github.com/JuliaDynamics/DynamicalBilliards.jl

# Features
The features of the DynamicalSystems.jl, as of version v1.6.1, are:

* Modular creation of a billiard table by putting together well-defined obstacles.
* Creation of random initial conditions in *any* user-created billiard table.
* Propagation of point particles in such billiard tables.
* Propagation of point particles in such billiard tables in magnetic fields.
* Exact calculation of the collision times, collision positions, velocities etc..
* Calculation of escape times of particles.
* Rich implementation of Ray-splitting billiards: A particle may propagate *through* an obstacle, given some arbitrary transmission probability.
* Calculation of Lyapunov exponents of trajectories (currently available only for propagation without magnetic fields).
* Flexible speed; users can trade between accuracy and speed.
* Extensive visualization library for plotting and animating particles and billiards.



# Description
Billiard systems have been used extensively in scientific research and played a
crucial role in the development of Chaos theory. A famous example is the Sinai billiard
(included in our package) which was one of the first low-dimensional systems proven to be ergodic [1].

Even though the study of billiard systems started decades ago, there are still new
surprises and plenty of research to be done with them. In [2] this is described perfectly, where Bunimovich and Vela-Arevalo (two pioneers in the field) summarize
new insights in the field of dynamical billiards.

Our package has plenty of applications, due to the large amount of features. It even offers the first time possibility
of implementing ray-splitting billiards, since, to the best of our knowledge, there is no other open sourced project that offers this. In addition,
the package is written in pure Julia [3] a new programming language with many advantages.
Because Julia is *dynamic*, it allows users to experiment with billiards flexibly and *interactively* and have access to a very easy-to-use high level interface.
On the other hand, because Julia is as fast as C, there is no worry about the speed of the algorithms.

Another advantage of using Julia is that the source code is clear, concise and very easy to understand even for beginners. This allows our package to be extendable in an almost trivial manner. In fact, in this [tutorial page](https://juliadynamics.github.io/DynamicalBilliards.jl/latest/tutorials/own_obstacle/) [4] we show users that all they need to do to define a completely new type of `Obstacle` (the things particles collide with) is only *four functions*, with an average of just 50 lines of code. Everything else is taken care of with the modular approach of our package and the power of abstraction that Julia provides.

DynamicalBilliards.jl calculates *exactly* all collisions between particles and obstacles, by solving 1st and 2nd order polynomial equations, which results to accuracies of the order of `1e-14` to `1e-16` (numerical limit for 64-bit floating point numbers). For this to work, particles are considered point particles. This advantage of exact solutions however also comes with a drawback: our package currently does not support (a) interactions between particles and (b) external electrostatic potentials.

We would like to bring into attention another new billiard package, called "Bill2D" [5], which is written in C instead, and supports interactions between particles. There is a some amount of overlap between DynamicalBilliards.jl and Bill2D, however both packages offer a lot of unique features not included in each other.


# References

[1] : Yakov G. Sinai, Russian Mathematical Surveys **25**, page 137 (1970), DOI: 10.1070/RM1970v025n02ABEH003794

[2] : L. Bunimovich & L. Vela-Arevalo, Chaos **25**, page 097614 (2015), DOI: 10.1063/1.4916330

[3] : https://julialang.org/, DOI: 10.1137/141000671

[4] : https://juliadynamics.github.io/DynamicalBilliards.jl/latest/tutorials/own_obstacle/

[5] : J. Solanpääa, P. J. J. Luukkob & E. Räsänena, Computer Physics Communications **199**, page 133 (2016), DOI: 10.1016/j.cpc.2015.10.006
