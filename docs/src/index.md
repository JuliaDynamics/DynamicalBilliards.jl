# Dynamical Billiards

`DynamicalBilliards.jl` is an easy-to-use, modular and extendable Julia package for Dynamical Billiard systems in two dimensions.

The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:
* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard. 
* **Ray-Splitting billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiard systems that have already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : functions for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more. Also includes animation related content.


## Installation

---

Since this package is not yet registered, use `Pkg.clone` until the registration is complete.

## Usage

---

For a crash course on how to use `DynamicalBilliards.jl`, you should visit the [Basic Usage](/basic/basic_usage) section.

If however, you want to make the most out if, the following tutorials offer detailed descriptions:
- [How to define your custom Billiard Table](/tutorials/billiard_table)
- [Using Ray-Splitting billiards](/tutorials/ray-splitting)
- [Creating your own Obstacle Type](/tutorials/own_obstacle)
- [Examples page](/tutorials/examples)

The [Library](/basic/library) section has the docstrings of all exported names in convenient groups.