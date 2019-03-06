# 3.5
* New function `visited_obstacles!`.
* 2D and 3D chaotic phase space volumes now in `MushroomTools`.

# 3.4
* Lyapunov computation and perturbation growth can now accept integer time (which means amount of collisions).

# 3.3
* Replaced all instances of `pertubation` with `perturbation`. A deprecation is thrown.

# 3.2
* New perturbation growth functions.

# 3.1
* `law_of_refraction` : new function that returns transmission and refraction functions based on Snell's law and refractive indices.
* Plenty of doc corrections.

# 3.0

## Breaking changes
* Plotting functions now overload `plot` instead. There is no more `plot_obstacle`, `plot_billiard` and `plot_particle`. `plot_boundarymap` and `animate_evolution` remain the same though.
* `construct` is removed.
* `boundarymap` call signature and return values are reworked.

## Enhancements / new features
* Brand new evolution animation function (replaced the old). It now animates in real-time!!! (you can use `dt = Inf` to obtain the old behavior). In addition it supports animating the evolution of multiple particles in parallel!!!
* New function `timeseries` creates the timeseries of the particle evolution in a billiard. This does pretty much what `construct(evolve...)` used to do.
* New obstacle: `Ellipse`.
* Much more robust propagation algorithm that is less prone to errors and "weird behaviors"!
* Test suite reworked almost from scratch: More tests, more specific tests, more robust tests, easier to debug tests!
* `totallength` is exported.
* `plot(bd, x, y)` (using the timeseries of `timeseries`) now also works for non-periodic billiards as well, for convenience.
* Automatic parallelization is now possible for some functions, given by the `parallelize` function.

## Low-Level changes
These changes are not actually breaking, unless someone used the low-level interface. The docs also changed and much less than the low level interface is exposed.

* Renamed `collisiontime` to just `collision`, since now it returns both the time and the estimated collision point.

* Many low-level functions are not exported any more, for safety and because it didn't make much sense: `normvalvec, distance, cellsize, ``propagate!`, `relocate!`, `resolvecollision!`, `periodicity!`, `specular!` `realangle`. For the low level interface only `propagate!`, `bounce!`, `collision` and `next_collision` are exposed. Only `bounce!` is exposed as public API. The low level interface is _still_ documented though.

* Change of the internal propagation algorithm:
  1. the function `collision` (previously `collisiontime`) returns both the time until collision *as well as* the collision point (most methods computed it already anyways).
  2. The `timeprec` and all those nonsense are removed. It is the duty of `collision` to return `Inf` if the particle is very close to the collision point (`realangle` uses `accuracy`, while for straight propagation the direction of travel is enough).
  3. The particle is then "teleported" to this point by setting its `pos`. For the case of magnetic propagation the velocity is propagated by the collision time. This brings very high accuracy, as the velocity vector is not multiplied with time and then added to the position vector.
  4. The `distance` is checked. If it has wrong sign, a `relocation` simply relocates the particle in the direction of the `normalvec`, for amount `(1+eps())*distance`.
  5. Specular reflection / periodicity is done.

* Renamed `distancecheck` to `accuracy`

# v2.3
* Now you can write `p.ω` as well as `p.omega` for magnetic particles.
* New `ispinned` function that returns Bool of whether a particle is pinned or not. Also works with periodic billiards.

# v2.2
* Standard billiards can now also be created with keyword arguments.
* Logo billiard is now exported by the function `billiard_logo`.

# v2.1
* Better limits of axis for periodic rectangular billiards.
* Added some methods for high level functions (like e.g. `evolve`) that if not given a particle they pick on from `randominside`.
* Documentation improvements.
* Bug fix in the expressions for the chaotic phase space volume of mushrooms.

# v2.0

## New Features!
* **3 orders of magnitude performance gains on all functions!!!**
  * Reduced a lot of allocations done all over the place. In most places allocations done are now exactly zero! ZEROOOOOOOOOOOO
  * Fixed many instances of broadcasting with static vectors (which is bad).
  * Utilized metaprogramming to manually unroll some loops.
  * If a billiard is periodic it is now known at compile time. This gives massive performance benefits in functions like `evolve`, which have to check if the collision occured with periodic walls.

* **Symver will now be properly respected from 2.0 onwards**.
* Automatic saving of animations to video via `ffmpeg`!!!
* Hexagonal periodic plotting.
* Lyapunov exponents for magnetic particles are now possible!
* Added boundary map computation function which works
  for any billiard and any particle. It assumes that the obstacles are
  sorted counter clockwise.
  * Added `to_bcoords`, `totallength`
  * Added `plot_boundarymap` that plots the boundary map and the obstacle boundaries.
* Added robust coordinate change transformation from 3D real space to 2D boundary space, see `to_bcoords`, `from_bcoords` and `arclengths`.
* Added two novel high-level functions that compute the phase-space portion an orbit covers as the particle evolves: `boundarymap_portion` and `phasespace_portion`.
* Added Poincare surface of section function, which computes intersections with arbitrary planes!
* It is now possible to affect many different obstacles during ray-splitting!
* Plotting is now available the moment the user does `using PyPlot`. Done through the `Requires` module. The function `enableplotting()` does not exist anymore!
* Re-organized all source code into a much more readable state, and as a result significantly reduced the total lines of code.
* added `evolve` function that simply deepcopies particle.
* Added convenience function to compute the mean collision time in a billiard.
* New low-level function `bounce!` that propagates a particle from one collision to the
  next. In essence does what `evolve!` does with `t=1`, but without creating a bunch
  of saving stuff. All high level functions use `bounce!`.

## Syntax and other breaking changes
* Default colors for plotting have been changed (random obstacles are purple,
  particle orbit is `"C0"`).
* **[BREAKING]** Overhauled what a "billiard table" is: Now, called simply `Billiard` is a
  dedicated struct that stores the obstacles as a tuple. This means that
  all functions do not accept a `Vector{Obstacle}` anymore but rather a `Billiard`.
* **[BREAKING]** `timeprec` now takes arguments `timeprec(::Particle, ::Obstacle)` to utilize better
  multiple dispatch and reduce code repetition.
* **[BREAKING]** `realangle` now only takes one intersection and simply returns the real angle.
* **[BREAKING]** `animate_evolution` now does not have `!` at the end, because it deepcopies the particle.
* **[BREAKING]** Re-worked ray-splitting: We now use the `RaySplitter` struct. See docs.

---

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
