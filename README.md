## DynamicalBilliards.jl
A Julia package for Dynamical billiard systems in two dimensions.
The goals of the package is to provide a flexible and intuitive framework for fast implementation of billiard systems of arbitrary construction.

The core of `DynamicalBilliards.jl` is separated in simple and cohesive modular structures:
* **Straight propagation** : The standard billiard dynamical system. A particle is propagating in a straight line, until a specular reflection is performed at a boundary.
* **Magnetic propagation** : Instead of a straight line, the orbit of the particle is a circle, like electrons in a perpendicular magnetic field. The particle still undergoes specular reflections at the boundaries of the billiard. 
* **Ray-splitting Billiards** : A semiclassical implementation of the dynamical billiard. After a collision of a particle with a boundary, the particle may propagate *through* the boundary given some arbitrary probability and transmission law.
* **Standard billiards** : A library of pre-constructed billiards that have been already been used in Physics/Mathematics (e.g. Sinai, periodic Sinai, Buminovich etc.)
* **Visualization** : routines for plotting and visualizing aspects of a billiard system, such as obstacles, orbits and more.
