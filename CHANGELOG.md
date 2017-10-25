# v1.6.1
Updated the documentation to reflect the new changes of v1.6.0

# v1.6.0 - WARNING: DOCUMENTATION NOT UPDATED. SEE v1.6.1
* **[BREAKING]** Function `animate_evolution` has been renamed to `animate_evolution!`
  to remind users that it mutates the particle.
* **[BREAKING]** `FiniteWall` has been renamed to `InfiniteWall` and works as before
  for convex billiards. A new type `FiniteWall` is introduced that can work for
  non-convex billiards.
  * `FiniteWall` has some extra fields for enabling this.
  * `FiniteWall` has a boolean field `isdoor`, that designates the given wall to be
    `Door`. This is used in `escapetime`.
* *MASSIVE*: Added function `escapetime(p, bt)` which calculates the escape time of a particle
  from a billiard table. The escape time is the time until the particle collides
  with a `Door` (any `Door`).
* `animate_evolution!` can create a new figure and plot the billiard table on
  user input. This happens by default.
* Bugfix: relocation in magnetic case was not adaptive (for the backwards method).
* *MASSIVE*: Added a `Semicircle` type! For both types of evolution!
    * added Bunimovich  billiard `billiard_bunimovich`
    * added mushroom billiard `billiard_mushroom`


# v1.5.0
* Added possibility to calculate the Lyapunov spectrum of a billiard
  system. Currently this is available only for `Particle`s.
    * use the function `lyapunovspectrum` for the computation.
* Changed the relocating algorithm to be geometric. i.e. the time adjusting is done
  geometrically (self-multiplying the adjustment by 10 at each repeated step).
* Magnetic propagation and straight propagation now both use `timeprec(T)` (see below).
  For both cases the maximum number of geometric iterations is 1 (for the default
  value of `timeprec(T)`.
* This `timeprec` cannot be used PeriodWall and RaySplitting obstacles with
  MagneticParticles because when magnetic and relocating forward you get extremely
  shallow angles and you need huge changes in time for even tiny changes in position.
  * For this special case the function `timeprec_forward(T)` is used instead. This
    function results to on average 3-5 geometric relocation steps.
* Fixed many issues and bugs.

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
