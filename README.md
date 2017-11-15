![DynamicalBilliards Logo: The Julia billiard](http://i.imgur.com/NKgzYrt.gif)

A Julia package for dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction.

| **Documentation**   | **Citation** | **Travis**     | **AppVeyor** | **Gitter** |
|:--------:|:--------:|:---------------:|:-----:|:-----:|
|[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://JuliaDynamics.github.io/DynamicalBilliards.jl/latest), [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/DynamicalBilliards.jl/stable)| [![status](http://joss.theoj.org/papers/753469f6b18c9c38127a7727d13c87cd/status.svg)](http://joss.theoj.org/papers/753469f6b18c9c38127a7727d13c87cd) | [![Build Status](https://travis-ci.org/JuliaDynamics/DynamicalBilliards.jl.svg?branch=master)](https://travis-ci.org/JuliaDynamics/DynamicalBilliards.jl) | [![Build status](https://ci.appveyor.com/api/projects/status/ksgb8pv5xl0j315y?svg=true)](https://ci.appveyor.com/project/JuliaDynamics/dynamicalbilliards-jl-wt09b) | [![Gitter](https://img.shields.io/gitter/room/nwjs/nw.js.svg)](https://gitter.im/JuliaDynamics/Lobby)


The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:

* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard.
* **Ray-splitting billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiard systems that have already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : Functions for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more. Also includes animation related content.
* **Lyapunov Spectrum** : Calculate the lyapunov spectrum of the trajectory of a particle in an arbitrary billiard table. Currently this is only available for `Particle`s.

**NOTICE:** This package does not support collision between particles (currently), since
all particles are considered point-particles.

## Installation
This package is registered, simply use `Pkg.add("DynamicalBilliards")` to install it.

The master branch of `DynamicalBilliards` is used for development purposes. It is not advised to use `Pkg.checkout("DynamicalBilliards")`, unless you want to contribute to the development of the package.

## Plotting
Plotting in `DynamicalBilliards` is done through the [`PyPlot` package]("https://github.com/JuliaPy/PyPlot.jl"). However, all plotting-related functions are not available by default but only "on-demand". Use `DynamicalBilliards.enableplotting()` to bring them into scope.

**WARNING**: You must be able to `using PyPlot` if you want to use the plotting capabilities of `DynamicalBilliards`! If you are having trouble installing `PyPlot` you can always use the minimal Python installation through miniconda by running these lines in your Julia terminal:

```julia
ENV["PYTHON"]=""; Pkg.add("PyCall"); Pkg.build("PyCall");
Pkg.add("PyPlot"); using PyPlot;
```

## Acknowledgements
This package is mainly developed by George Datseris. However, this development would not have been possible without significant help from other people:

1. [Diego Tapias](https://github.com/dapias) (@dapias) Contributed the lyapunov spectrum calculation methods.
1. [David. P. Sanders](https://github.com/dpsanders) (@dpsanders) and [Ragnar Fleischmann](https://www.ds.mpg.de/person/20199/118124) contributed in fruitful discussions about the programming and physics of Billiard systems all-around.
2. [Christopher Rackauckas](https://github.com/ChrisRackauckas) (@ChrisRackauckas) helped set-up the continuous integration, testing, documentation publishing and all around package development-related concepts.
3. [Tony Kelman](https://github.com/tkelman) (@tkelman) helped significantly in the package publication process, especially in making it work correctly without destroying METADATA.jl.
