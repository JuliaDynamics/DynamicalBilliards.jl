
Below you find all the docstrings of all exported names of `DynamicalBilliards.jl` in convenient groups.


<a id='Particles-1'></a>

## Particles

<a id='DynamicalBilliards.Particle' href='#DynamicalBilliards.Particle'>#</a>
**`DynamicalBilliards.Particle`** &mdash; *Type*.



```
Particle <: AbstractParticle
```

Two-dimensional particle in a billiard table.

**Fields:**

  * `pos::SVector{2,Float64}` : Current position vector.
  * `vel::SVector{2,Float64}` : Current velocity vector (always of measure 1).
  * `current_cell::SVector{2,Float64}` : Current "cell" the particle is located at. (Used only in periodic billiards)

**Additional constructors:**

```julia
Particle{T<:Real}(ic::Vector{T}) #where ic = [x0, y0, φ0]
Particle(x::Real, y::Real, φ::Real)
Particle() = Particle(rand(), rand(), rand()*2π)
```

<a id='DynamicalBilliards.MagneticParticle' href='#DynamicalBilliards.MagneticParticle'>#</a>
**`DynamicalBilliards.MagneticParticle`** &mdash; *Type*.



```
MagneticParticle <: AbstractParticle
```

Two-dimensional particle in a billiard table with perpendicular magnetic field.

**Fields:**

  * `pos::SVector{2,Float64}` : Current position vector.
  * `vel::SVector{2,Float64}` : Current velocity vector (always of measure 1).
  * `current_cell::SVector{2,Float64}` : Current "cell" the particle is located at. (Used only in periodic billiards)
  * `omega::Float64` : Angular velocity (cyclic frequency) of rotational motion. Radius of rotation is `r=1/omega`.

**Additional constructors:**

```julia
MagneticParticle{T<:Real}(ic::Vector{T}, ω::Real) #where ic = [x0, y0, φ0]
MagneticParticle(x0::Real, y0::Real, φ0::Real, ω::Real)
MagneticParticle() = MagneticParticle([rand(), rand(), rand()*2π], 1.0)
```

<a id='DynamicalBilliards.magnetic2standard' href='#DynamicalBilliards.magnetic2standard'>#</a>
**`DynamicalBilliards.magnetic2standard`** &mdash; *Function*.



```julia
magnetic2standard(p::MagneticParticle, use_cell = true)
```

Create a standard `Particle` from a `MagneticParticle`.

<a id='DynamicalBilliards.standard2magnetic' href='#DynamicalBilliards.standard2magnetic'>#</a>
**`DynamicalBilliards.standard2magnetic`** &mdash; *Function*.



```julia
standard2magnetic(omega, p::Particle, use_cell = true)
```

Create a `MagneticParticle` from a `Particle`.

<a id='DynamicalBilliards.cyclotron' href='#DynamicalBilliards.cyclotron'>#</a>
**`DynamicalBilliards.cyclotron`** &mdash; *Function*.



```julia
cyclotron(p::MagneticParticle, use_cell = false)
```

Return center and radius of circular motion performed by the particle based on `p.pos` (or `p.pos + p.current_cell`) and `p.vel`.


<a id='Obstacles-1'></a>

## Obstacles

<a id='DynamicalBilliards.Obstacle' href='#DynamicalBilliards.Obstacle'>#</a>
**`DynamicalBilliards.Obstacle`** &mdash; *Type*.



```
Obstacle
```

Obstacle supertype.

<a id='DynamicalBilliards.Disk' href='#DynamicalBilliards.Disk'>#</a>
**`DynamicalBilliards.Disk`** &mdash; *Type*.



```
Disk <: Circular
```

Disk-like obstacle with propagation allowed outside of the circle.

**Fields:**

  * `c::SVector{2,Float64}` : Center.
  * `r::Float64` : Radius.
  * `name::String` : Some name given for user convenience.

Constructors accept any vectors convertible to SVector{2,Float64}.

<a id='DynamicalBilliards.Antidot' href='#DynamicalBilliards.Antidot'>#</a>
**`DynamicalBilliards.Antidot`** &mdash; *Type*.



```
Antidot <: Circular
```

Disk-like obstacle that allows propagation both inside and outside of the disk. Used in ray-splitting billiards.

**Fields:**

  * `c::SVector{2,Float64}` : Center.
  * `r::Float64` : Radius.
  * `where::Bool` : Flag that keeps track of where the particle is currently propagating. `true` stands for *outside* the disk, `false` for *inside* the disk.
  * `name::String` : Name of the obstacle given for user convenience.

Constructors accept any vectors convertible to SVector{2,Float64}.

<a id='DynamicalBilliards.FiniteWall' href='#DynamicalBilliards.FiniteWall'>#</a>
**`DynamicalBilliards.FiniteWall`** &mdash; *Type*.



```
FiniteWall <: Wall
```

Wall obstacle imposing specular reflection during collision.

**Fields:**

  * `sp::SVector{2,Float64}` : Starting point of the Wall.
  * `ep::SVector{2,Float64}` : Ending point of the Wall.
  * `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the particle *will come from before a collision* (pointing towards the inside the billiard table). The size of the vector is irrelevant.
  * `name::String` : Name of the obstacle, e.g. "left wall", given for user convenience.

Constructors accept any vectors convertible to SVector{2,Float64}.

<a id='DynamicalBilliards.PeriodicWall' href='#DynamicalBilliards.PeriodicWall'>#</a>
**`DynamicalBilliards.PeriodicWall`** &mdash; *Type*.



```
PeriodicWall <: Wall
```

Wall obstacle that imposes periodic boundary conditions upon collision.

**Fields:**

  * `sp::SVector{2,Float64}` : Starting point of the Wall.
  * `ep::SVector{2,Float64}` : Ending point of the Wall.
  * `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the particle *will come from* (to the inside the billiard table). The size of the vector is **important**. This vector is added to a particle's `pos` during collision. Therefore the size of the normal vector must be correctly associated with the size of the periodic cell.
  * `name::String` : Name of the obstacle, e.g. "left boundary", given for user convenience.

Constructors accept any vectors convertible to SVector{2,Float64}.

<a id='DynamicalBilliards.SplitterWall' href='#DynamicalBilliards.SplitterWall'>#</a>
**`DynamicalBilliards.SplitterWall`** &mdash; *Type*.



```
SplitterWall <: Wall
```

Wall obstacle imposing specular reflection during collision.

**Fields:**

  * `sp::SVector{2,Float64}` : Starting point of the Wall.
  * `ep::SVector{2,Float64}` : Ending point of the Wall.
  * `normal::SVector{2,Float64}` : Normal vector to the wall, pointing to where the particle *will come from before a collision*. The size of the vector is irrelevant.
  * `where::Bool` : Flag that keeps track of where the particle is currently propagating. `true` is associated with the `normal` vector the wall is instantiated with.
  * `name::String` : Name of the obstacle, e.g. "ray-splitting wall 1", given for user convenience.

Constructors accept any vectors convertible to SVector{2,Float64}.

<a id='DynamicalBilliards.normalvec' href='#DynamicalBilliards.normalvec'>#</a>
**`DynamicalBilliards.normalvec`** &mdash; *Function*.



```julia
normalvec(obst::Obstacle, position)
```

Return the vector normal to the obstacle at the given position (which is assumed to be very close to the obstacle's boundary).

The normal vector of any Obstacle must be looking towards the direction a particle is expected to come from.

<a id='DynamicalBilliards.distance' href='#DynamicalBilliards.distance'>#</a>
**`DynamicalBilliards.distance`** &mdash; *Function*.



```julia
distance(p::AbstractParticle, o::Obstacle)
```

Return the **signed** distance between particle `p` and obstacle `o`, based on `p.pos`. Positive distance corresponds to the particle being on the *allowed* region of the Obstacle. E.g. for a `Disk`, the distance is positive when the particle is outside of the disk, negative otherwise.

```julia
distance(p::AbstractParticle, bt::Vector{Obstacle})
```

Return minimum `distance(p, obst)` for all `obst` in `bt`, which can be negative.

<a id='DynamicalBilliards.randominside' href='#DynamicalBilliards.randominside'>#</a>
**`DynamicalBilliards.randominside`** &mdash; *Function*.



```julia
randominside(bt::Vector{Obstacle}[, omega])
```

Return a particle with correct (allowed) initial conditions inside the given billiard table defined by the vector `bt`. If supplied with a second argument the type of the returned particle is `MagneticParticle`, with angular velocity `omega`. Else, it is `Particle`.


<a id='Propagation-1'></a>

## Propagation

<a id='DynamicalBilliards.resolvecollision!' href='#DynamicalBilliards.resolvecollision!'>#</a>
**`DynamicalBilliards.resolvecollision!`** &mdash; *Function*.



```
resolvecollision!(p::AbstractParticle, o::Obstacle)
```

Resolve the collision between particle `p` and obstacle `o`. If the obstacle is not a periodic wall, the function performs specular reflection. If it is a periodic wall, it performs the periodicity condition.

`resolvecollision!()` takes special care so that the particle is always inside the correct side of the billiard table, in order to avoid particle leakage. Specifically, it calculates the distance from particle and obstacle and, depending on the obstacle type, makes necessary adjustments by propagating the particle forwards or backwards in time using **linear** motion.

```
resolvecollision!(p, o, T::Function, θ::Function, new_ω::Function)
```

This is the ray-splitting implementation. The three functions given are drawn from the ray-splitting dictionary that is passed directly to `evolve!()`. For a calculated incidence angle φ, if T(φ) > rand(), ray-splitting occurs. (See the section "Ray-Splitting" of the official documentation for more info.)

<a id='DynamicalBilliards.collisiontime' href='#DynamicalBilliards.collisiontime'>#</a>
**`DynamicalBilliards.collisiontime`** &mdash; *Function*.



```julia
collisiontime(p::AbstractParticle, o::Obstacle)
```

Calculate the collision time between given particle and obstacle.

The funtion chooses the appropriate method depending on the type of particle (magnetic or not) as well as the type of the obstacle. Returns the time that the particle, given its current position and Type, must be propagated to reach the collision point. This time can be given directly to `propagate!(p, time)` which brings the particle to the collision point.

In the case of magnetic propagation, there are always two possible collisions. The function internally decides which of the two will occur first, based on the sign of the angular velocity of the magnetic particle.

<a id='DynamicalBilliards.propagate!' href='#DynamicalBilliards.propagate!'>#</a>
**`DynamicalBilliards.propagate!`** &mdash; *Function*.



```julia
propagate!(p::AbstractParticle, t)
```

Propagate the particle `p` for given time `t`, changing appropriately the the `p.pos` and `p.vel` fields.

For a `Particle` the propagation is a straight line (i.e. velocity vector is constant).

For a `MagneticParticle` the propagation is circular motion with cyclic frequency `p.omega` and radius `1/p.omega`.

<a id='DynamicalBilliards.evolve!' href='#DynamicalBilliards.evolve!'>#</a>
**`DynamicalBilliards.evolve!`** &mdash; *Function*.



```
evolve!(p::AbstractParticle, bt::Vector{Obstacle}, ttotal)
```

Evolve the given particle `p` inside the billiard table `bt` for a total amount of time `ttotal`. Return the states of the particle between collisions.

The evolution takes into account the particle's Type. E.g. if `typeof(p) == MagneticParticle` then magnetic evolution will take place.

**Calling**

Call the function like:

```julia
t, poss, vels, (ω)* = evolve!(p, bt, ttotal)
```

(see "Returns" section for more)

To get the position, velocity and time timeseries from the above output, use the function `construct`: `xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ttotal)...)`

**Returns**

As noted by the "!" at the end of the function, it changes its argument `p` (particle). Most importantly however, this function also returns the main output expected by a billiard system. This output is a tuple of 3 (or 4) vectors:

  * `ct::Vector{Float64}` : Collision times.
  * `poss::Vector{SVector{2}}` : Positions during collisions.
  * `vels:: Vector{SVector{2}})` : Velocities **exactly after** the collisions.
  * `ω`, either `Float64` or `Vector{Float64}` : Angular velocity(/ies).

In the case of straight propagation, only the first 3 are returned.

In the case of magnetic propagation, the 4th value is returned as well. This is either the angular velocity of the particle (`Float64`), or in the case of ray-splitting it is a vector of the angular velocities at each time step (`Vector`).

The time `ct[i]` is the time necessary to reach state `poss[i+1], vels[i+1]` starting from the state `poss[i], vels[i]`. That is why `ct[1]` is always 0 since `poss[1], vels[1]` are the initial conditions. The angular velocity `ω[i]` is the one the particle has while propagating from state `poss[i], vels[i]` to `i+1`.

Notice that at any point, the velocity vector `vels[i]` is the one obtained **after** the specular reflection of the (i-1)th collision. The function `construct` takes that into account.

**Ray-splitting billiards**

No matter how complex ray-splitting processes you want, and irrespectively of how many obstacles in the billiard table can perform ray-splitting, there is only a single difference on the main function call: The `evolve!()` function is supplemented with a fourth argument, `ray_splitter::Dict{Int, Vector{Function}}`. This dictionary object handles all ray-splitting processes in the billiard system. It is a map of the Obstacle index within the billiard table to the ray-splitting functions: (φ is the angle of incidence)

  * T(φ, where, ω) : Transmission probability.
  * θ(φ, where, ω) : Transmission (aka refraction) angle.
  * ω_new(ω, where) : Angular velocity after transmission.

For more information and instructions on defining these functions please visit the official documentation.

<a id='DynamicalBilliards.construct' href='#DynamicalBilliards.construct'>#</a>
**`DynamicalBilliards.construct`** &mdash; *Function*.



```julia
construct(ct, poss, vels[, ω][, dt=0.01])
```

Given the main output of this package (see `evolve!()` function) construct the timeseries of the position and velocity, as well as the time vector.

In case of not given ω (or ω == 0), straight construction takes place. In case of ω != 0 or in case of ω::Vector{Real} magnetic construction takes place.

The additional optional argument of `dt` (only valid for Magnetic construction) is the timestep with which the timeseries are constructed.

**Calling**

Call this function like:

```julia
# Straight propagation:
xt, yt, vxt, vyt, ts = construct(ct, poss, vels)
# Magnetic propagation:
xt, yt, vxt, vyt, ts = construct(ct, poss, vels, ω, dt)
# Any-kind-of propagation (skip dt for straight):
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, ttotal[, ray_splitter])..., dt)
```

There is no difference in the `construct` call for Ray-Splitting billiards if one uses the awesome ellipsis operator of Julia. The reason for that is that `evolve!()` also returns the vector of angular velocities when necessary.

**Returns**

A tuple of the following:

  * xt::Vector{Float64} : x position time-series
  * yt::Vector{Float64} : y position time-series
  * vxt::Vector{Float64} : x velocity time-series
  * vyt::Vector{Float64} : y velocity time-series
  * ts::Vector{Float64} : Continuous time vector

<a id='DynamicalBilliards.isphysical' href='#DynamicalBilliards.isphysical'>#</a>
**`DynamicalBilliards.isphysical`** &mdash; *Function*.



```
isphysical(raysplitter::Dict{Int, Vector{Function}}; only_mandatory = false)
```

Return `true` if the given ray-splitting dictionary represends the physical world.

Specifically, check if (φ is the incidence angle):

  * Critical angle means total reflection: If θ(φ) ≥ π/2 then T(φ) = 0
  * Transmission probability is even function: T(φ) ≈ T(-φ)
  * Refraction angle is odd function: θ(φ) ≈ -θ(-φ)
  * Ray reversal is true: θ(θ(φ, where, ω), !where, ω) ≈ φ
  * Magnetic conservation is true: (ω_new(ω_new(ω, where), !where) ≈ ω

The first property is mandatory and must hold for correct propagation. They keyword `only_mandatory` notes whether the rest of the properties should be tested or not.


<a id='Standard-Billiards-1'></a>

## Standard Billiards

<a id='DynamicalBilliards.billiard_rectangle' href='#DynamicalBilliards.billiard_rectangle'>#</a>
**`DynamicalBilliards.billiard_rectangle`** &mdash; *Function*.



```
billiard_rectangle(x=1.0, y=1.0; periodic = false)
```

Return a vector of obstacles that defines a rectangle billiard of size (`x`, `y`).

<a id='DynamicalBilliards.billiard_sinai' href='#DynamicalBilliards.billiard_sinai'>#</a>
**`DynamicalBilliards.billiard_sinai`** &mdash; *Function*.



```
billiard_sinai(r, x=1.0, y=1.0; periodic = false)
```

Return a vector of obstacles that defines a Sinai billiard of size (`x`, `y`) with a disk in its center, of radius `r`.

In the periodic case, the system is also known as "Lorentz Gas".

<a id='DynamicalBilliards.billiard_polygon' href='#DynamicalBilliards.billiard_polygon'>#</a>
**`DynamicalBilliards.billiard_polygon`** &mdash; *Function*.



```
billiard_polygon(n::Int, R, center = [0,0]; periodic = true)
```

Return a vector of obstacles that defines a regular-polygonal billiard table with `n` sides, radius `r` and given `center`. If `n` is even, you may choose a periodic version of the billiard.

Note: `R` denotes the so-called outer radius, not the inner one.

