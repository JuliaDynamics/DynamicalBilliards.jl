![DynamicalBilliards Logo: The Julia billiard](http://i.imgur.com/NKgzYrt.gif)

`DynamicalBilliards` is an easy-to-use, modular and extendable Julia package for
dynamical billiards in two dimensions.

!!! warn "Version 2.0 incoming!"
    We are currently developing 2.0 of `DynamicalBilliards` which will bring massive
    performance increase, a lot of new features and unfortuantelly some breaking changes! You can view the [changelog](https://github.com/JuliaDynamics/DynamicalBilliards.jl/blob/master/CHANGELOG.md) to prepare for up coming changes.


!!! tip "Julia Billiard animation"
    Check out the example in the [tutorials](tutorials/examples/#julia-logo-billiard) page to see the code that created and animated the "Julia Billiard"!

## Introduction

A dynamical billiard is a system where a particle is propagating as a straight line from obstacle to obstacle, performing specular reflection at the boundary of the obstacles. Billiard systems have been used extensively in mathematics, nonlinear dynamics and chaos and played an important role in the development of nonlinear science.

The [wikipedia page](https://en.wikipedia.org/wiki/Dynamical_billiards) has many examples of different types of billiards. Also, the [scholarpedia](http://www.scholarpedia.org/article/Dynamical_billiards) entry is a good read on the subject.

## Features

* Modular creation of a billiard from well defined obstacles
* Straight propagation of a particle in a billiard table
* Support for creating random initial conditions in an arbitrary
  billiard table
* Magnetic propagation, where the particle moves in a circle instead
  of a straight line (works with *any* billiard)
* Ray-splitting implementation: a particle may propagate
  through an obstacle given arbitrary transmission and refraction
  laws. This is also known as a "semiclassical billiard"
* Computation of Poincaré surfaces of section (also known as boundary maps) for any table and any particle
* Escape times
* Easy to use low-level interface
* Full support for visualizing and animating billiards and motion in billiards
* Brutal tests that confirm the package works and overcomes numerical precision issues

This package does not support finite-sized particles and, as a result, there is
also no support for collision between particles.

---


## Usage
It is highly suggested to first read the [High Level API](/basic/high_level).

See the [Visualizing](visualizing) page for plotting, animating, etc.

For more advanced usage see the [Low Level API](/basic/low_level).

The following tutorials offer detailed descriptions for various aspects of `DynamicalBilliards`:

- [How to define your custom Billiard](/tutorials/billiard_table)
- [Using Ray-Splitting billiards](/tutorials/ray-splitting)
- [Creating your own Obstacle Type](/tutorials/own_obstacle)
- [Examples page](/tutorials/examples)

The discussion about [numerical precision](/physics/#numerical-precision) is done in the [Physics](/physics) page.

The [Lyapunov Exponents](/lyapunovs) page has info on how to compute the Lyapunov
spectrum of a billiard.


---

## Citing
If you have used this package for research that resulted in a publication, please be
kind enough to cite the paper associated with DynamicalBilliards.jl. The DOI is
https://doi.org/10.21105/joss.00458 and you can cite as:

>G. Datseris, [The Journal of Open Source Software **2**, 458
(2017)](https://doi.org/10.21105/joss.00458).

or if you use BibTeX:
```
@article{Datseris2017,
  doi = {10.21105/joss.00458},
  url = {https://doi.org/10.21105/joss.00458},
  year  = {2017},
  month = {nov},
  volume = {2},
  number = {19},
  pages = {458},
  author = {George Datseris},
  title = {{DynamicalBilliards}.jl: An easy-to-use,  modular and extendable Julia package for Dynamical Billiard systems in two dimensions.},
  journal = {The Journal of Open Source Software}
}
```

## Installation

This package is registered, simply use `Pkg.add("DynamicalBilliards")` to install it.
The [stable documentation](https://juliadynamics.github.io/DynamicalBilliards.jl/stable/) accompanies the version installed with `Pkg.add()`.

Plotting is done through the `PyPlot` module. All plotting functions are brought
into scope when `using PyPlot` is done.
---

---
## Support
If you are having problems with `DynamicalBilliards.jl` do not hesitate to seek for support! There are numerous ways to do that:

1. Visit our [official chatroom](https://gitter.im/JuliaDynamics/Lobby) on Gitter: https://gitter.im/JuliaDynamics/Lobby
2. Open a new issue at our [GitHub issues page](https://github.com/JuliaDynamics/DynamicalBilliards.jl/issues).

---
## Contributing
Everyone is welcomed to contribute to `DynamicalBilliards.jl`! If you have some new
algorithm, types of Obstacles or anything new to add, do not hesitate! For formal
questions about e.g. structuring of code it is best to contact us through the [gitter
chatroom](https://gitter.im/JuliaDynamics/Lobby) or by opening a new Pull Request and asking for a review of your code.

If you would like to help but do not have anything new to contribute, please go ahead
and take a look at the [GitHub issues page](https://github.com/JuliaDynamics/DynamicalBilliards.jl/issues) of the package.
Some of the existing issues are easy to solve and are there specifically for people that would
like to contribute.

---
