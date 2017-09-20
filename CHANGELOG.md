# v1.5.0
* Added possibility to calculate the Lyapunov spectrum of a billiard
  system. Currently this is available only for `Particle`s.
    * use the function `lyapunovspectrum` for the computation.

# v1.4.0
* All major types defined by `DynamicalBilliards` have been changed to
  parametric types with parameter `T<:AbstractFloat`.
    * The above makes the billiard table type annotation be
      `bt::Vector{<:Obstacle{T}}) where {T<:AbstractFloat}` instead of
      the old `bt::Vector{Obstacle}`.
* Particle evolution algorithm has fundamentally changed.
  The way the algorithm works now is described in the documentation in the [Physics](https://juliadynamics.github.io/DynamicalBilliards.jl/latest/physics/#numerical-precision) page.
    * This point with conjuction with the above made the package **much faster**.
* Positional argument `warning` of `evolve!()` has been changed to **keyword argument**.
* The raysplitting functions must always accept 3 arguments, even in the case
  of straight propagation.
  The best way is to have the third argument have a default value.
* All `distance` functions can now take a position as an argument (giving a particle
  simply passes the position).
* Removed redundant functions like `magnetic2standard`.
* The package is now precompiled.
* Tests have been restructured to be faster, more efficient, and use the
  `@testset` type of syntax.
* Documentation has been made more clear.

# v1.3.0
Changelog will be kept from this version. See the [releases](https://github.com/JuliaDynamics/DynamicalBilliards.jl/releases) page for info on previous versions.
