![DynamicalBilliards v3.0 Logo: The Julia billiard](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/billiards/DynamicalBilliards_logo_animated.gif?raw=true)

`DynamicalBilliards` is an easy-to-use, modular and extendable Julia package for
dynamical billiards in two dimensions.
It is part of [JuliaDynamics](https://juliadynamics.github.io/JuliaDynamics/), an organization dedicated to creating high quality scientific software.

!!! info "Latest news"
    The [`timeseries!`](@ref) function now supports evolving particles until a certain condition (specified by a function) is met.

## About Billiards

A dynamical billiard is a system where a particle is propagating inside a domain,
bouncing from obstacle to obstacle (i.e. the *boundary* of the domain) by a specular reflection at the boundary of the obstacles. This basic idea can be extended in many ways, one of which is replacing the particle orbit from a straight line to a circle.

Billiard systems have been used extensively in mathematics, nonlinear dynamics and chaos and played an important role in the development of nonlinear science.
The [wikipedia page](https://en.wikipedia.org/wiki/Dynamical_billiards) has many examples of different types of billiards. Also, the [scholarpedia](http://www.scholarpedia.org/article/Dynamical_billiards) entry is a good read on the subject.

## Features

* Particles are evolved by solving *exactly* the geometric equations for intersections between lines, circles, ellipses, and other shapes. There are no approximations done regarding the dynamics.
* Modular creation of a [Billiard](@ref) from well defined obstacles. Arbitrary billiard shapes can be made and no shape is "hard coded".
* Full support for both *straight*  and *magnetic* propagation of a particle in a billiard table.
  * During magnetic propagation the particle orbit is a circle instead of a line!
  * All features exist for both types of propagation!
  * See the [High Level API](@ref) to get started!

* Support for creating [Random initial conditions](@ref) in an arbitrary billiard.
* [Ray-Splitting](@ref): a particle may propagate
  through an obstacle given arbitrary transmission and refraction
  laws. This is also known as a "semiclassical billiard".
* [Poincaré Sections](@ref) (intersections with arbitrary plane).
* [Boundary Maps](@ref).
* [Escape Times](@ref) & [Mean Collision Times](@ref).
* [Lyapunov Exponents](@ref).
* Support for both coordinate systems: 3D real space and boundary coordinates.
* Novel algorithms that compute the portion of either the 2D boundary space or the 3D real space that an orbit covers as a particle evolves. See the [Phase Spaces](@ref) section.
* Easy to use low-level interface, described at the [Internals](@ref) page.
* Specialized tools for mushroom billiards.
* Full support for [Visualizations and Animations for Billiards](@ref).
* Brutal tests that confirm the package works and overcomes numerical precision issues.

This package does not support finite-sized particles and, as a result, there is also no support for collision between particles.

## Citing
If you have used this package for research that resulted in a publication, please be
kind enough to cite the software paper associated with `DynamicalBilliards`.
The DOI is
https://doi.org/10.21105/joss.00458 and you can cite as:

> G. Datseris, [The Journal of Open Source Software **2**, 458 (2017)](https://doi.org/10.21105/joss.00458).

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

In addition, if you are using the functionality to compute Lyapunov exponents
in billiards, then please cite the following Chaos publication:
```
@article{Datseris2019,
  doi = {10.1063/1.5099446},
  url = {https://doi.org/10.1063/1.5099446},
  year = {2019},
  month = sep,
  publisher = {{AIP} Publishing},
  volume = {29},
  number = {9},
  pages = {093115},
  author = {George Datseris and Lukas Hupe and Ragnar Fleischmann},
  title = {Estimating Lyapunov exponents in billiards},
  journal = {Chaos: An Interdisciplinary Journal of Nonlinear Science}
}
```



## Installation

This package is registered, simply use `]` to get into the package manager mode and then type `add DynamicalBilliards` to install it.
The [stable documentation](https://juliadynamics.github.io/DynamicalBilliards.jl/stable/) accompanies the version installed with `add`.

To confirm the validity of your installation you can run the tests of `DynamicalBilliards`. This can be done via `] test DynamicalBilliards`.


## How to easily code a Billiard

We have created a [Jupyter notebook](https://nbviewer.jupyter.org/github/JuliaDynamics/JuliaDynamics/blob/master/tutorials/Billiards%20Example/billiards_example.ipynb) that showcases how easy it is to simulate a dynamical billiard using Julia. This [notebook](https://nbviewer.jupyter.org/github/JuliaDynamics/JuliaDynamics/blob/master/tutorials/Billiards%20Example/billiards_example.ipynb) is an educative example of both using Multiple Dispatch and of how the internal code of `DynamicalBilliards` works. It also highlights the extendibility of the core code.

Besides a Jupyter notebook, this tutorial is also available as an interactive article in the new platform NextJournal [here](https://nextjournal.com/julia/billiard).




## Support
If you are having any kind of problems with `DynamicalBilliards` do not hesitate to seek for support! There are numerous ways to do that:

1. Visit our [official chatroom](https://gitter.im/JuliaDynamics/Lobby) on Gitter: https://gitter.im/JuliaDynamics/Lobby
2. Open a new issue at our [GitHub issues page](https://github.com/JuliaDynamics/DynamicalBilliards.jl/issues).


## Contributing
Everyone is welcomed to contribute to `DynamicalBilliards`! If you have some new
algorithm, types of Obstacles or anything new to add, do not hesitate! For formal
questions about e.g. structuring of code it is best to contact us through the [gitter
chatroom](https://gitter.im/JuliaDynamics/Lobby) or by opening a new Pull Request and asking for a review of your code.

If you would like to help but do not have anything new to contribute, please go ahead
and take a look at the [GitHub issues page](https://github.com/JuliaDynamics/DynamicalBilliards.jl/issues) of the package.
Some of the existing issues are easy to solve and are there specifically for people that would
like to contribute.
