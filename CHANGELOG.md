# v2.0.0
* Added possibility to calculate the Lyapunov spectrum of a billiard
  system. Currently this is available only for `Particle`s.
* All major types defined by `DynamicalBilliards` have been changed to
  parametric types with parameter `T<:AbstractFloat`.
  * The above makes the billiard table type annotation be
    `bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}` instead of
    the old `bt::Vector{Obstacle}`.
* Positional argument `warning` of `evolve!()` has been changed to **keyword argument**.
* The raysplitting functions must always accept 3 arguments, even in the case
  of straight propagation. The best way is to hagve the third argument have a default value.
* All `distance` functions can now take a position as an argument (giving a particle
  simply passes the position).
* Removed redundant functions like `magnetic2standard`.
* The package is now precompiled.

* Here I talk about the new way propagation is handled. `relocate!` is super changed.

##not done yet:
* Tests have been restructured to be faster, more efficient, and use the
  `@testset` type of syntax.

# v1.3.0
Changelog will be kept from this version. See the [releases](https://github.com/JuliaDynamics/DynamicalBilliards.jl/releases) page for info on previous versions.
