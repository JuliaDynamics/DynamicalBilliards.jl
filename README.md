# DynamicalBilliards.jl

[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://Datseris.github.io/DynamicalBilliards.jl/stable)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://Datseris.github.io/DynamicalBilliards.jl/latest)
[![Build Status](https://travis-ci.org/Datseris/DynamicalBilliards.jl.svg?branch=master)](https://travis-ci.org/Datseris/DynamicalBilliards.jl)

A Julia package for dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction. ![Example animation](https://github.com/Datseris/DynamicalBilliards.jl/blob/master/images/plot_example.gif "Evolution of particle in a magnetic field.")

*(the code that generated this animation is shown at the end of this README)*

The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:
* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard.
* **Ray-splitting billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiard systems that have already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : functions for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more. Also includes animation related content.

NOTICE: This package does not support collision between particles. All particles are considered point-particles for all simulations offered by `DynamicalBilliards.jl`.

Documentation: coming soon.

*(all exported names of DynamicalBilliards.jl have detailed docstrings. Use `?` when in doubt)*

### Installation
This package is not yet registered. Use `Pkg.clone("https://github.com/Datseris/DynamicalBilliards.jl")` in order to install it.
