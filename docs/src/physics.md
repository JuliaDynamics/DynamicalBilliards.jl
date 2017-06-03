This page briefly discusses physical aspects of billiard systems.

## Implementation
Firstly one defines a billiard table and (if desired) the [ray-splitting dictionary](/tutorials/ray-splitting/#ray-splitter-dictionary). Then one creates a particle inside the defined billiard table. The algorithm followed for the propagation of a particle is the following:

1. Calculate the `collisiontime()` for **all** obstacles in the billiard table.
2. Choose the smallest time and `propagate!()` the particle for that amount of time. The obstacle corresponding to that minimum time is `colobst`.
3. `resolvecollision!()` between the particle and `colobst`:
    1. Check whether there is transmission or not (only for ray-splitting): `T(φ) > rand()`
    1. `relocate!()` the particle accordingly so that it is on the correct side of the billiard table.
    2. For no transmission (or no ray-splitting), perform `specular!()` reflection or `periodicity!()` conditions.
    3. Otherwise, implement the ray-splitting algorithm (not discussed here).  

5. Continue the loop 1-3 for a given amount of time.

In the standard billiard case, one can always exclude the collision with the previous obstacle. However, in both magnetic or ray-splitting cases this is not true anymore. Therefore the same algorithm is applied on all 3 cases for the sake of simplicity.

Notice that the `relocate!()` step is actually very important because it takes care that there is not particle "leakage": particles being outside the billiard table due to the finite precision of floating numbers.

!!! danger "Inefficiency "
    The first step of the algorithm is **very inefficient** for large billiard tables, since many of the obstacles could be safely excluded from search. If you have better solutions on how to find the minimum collision time you are more than welcome to submit a Pull Request or open an issue to discuss about it. Unfortunately, in the present literature,
    such an algorithm exists only for the periodic Sinai billiard specifically.

!!! note "Speed vs. Accuracy"
    `DynamicalBilliards.jl` was made under the consideration of as much accuracy as possible.
    There are numerous actions performed internally such that
    the particle propagation is accurate (bounded by `Float64` precision of course), and will not result in any kind of unexpected behavior,
    like for example particle leakage.
    Of course, this is a direct trade-off in the expense of speed.

## Ray-Splitting Functions
If `T` is the transmission probability function, then the condition for transmission is simply: `T(φ, pflag, ω) > rand()`. If it returns `true`, transmission (i.e. ray-splitting) will happen. As it has already been discussed in the [Ray-Splitting tutorial](/tutorials/ray-splitting), the condition of total internal reflection must be taken care of by the user.

The [three key functions](/tutorials/ray-splitting/#ray-splitting-functions) given to the `ray-splitter` dictionary must have some properties in order to have physical meaning, like for example that the scattering probability function is even towards φ. One of these properties is absolutely **mandatory** for this package to work properly. This is the property of total internal reflection, i.e. if the refraction angle is calculated to be greater/equal than π/2, no transmission can happen. **This condition is not assured internally** and therefore you must be sure that your transmission probability function satisfies it. In the above example, the function `T` makes sure to return 0 in that case.

In order to test if the `raysplitter` dictionary you have defined has physical meaning, the function `isphysical()` is provided:

```julia
    isphysical(raysplitter::Dict{Int, Any}; only_mandatory = false)
```
Return `true` if the given ray-splitting dictionary has physically plausible properties.

Specifically, check if (φ is the incidence angle, θ the refraction angle):
* Critical angle means total reflection: If θ(φ) ≥ π/2 then T(φ) = 0
* Transmission probability is even function: T(φ) ≈ T(-φ) at ω = 0
* Refraction angle is odd function: θ(φ) ≈ -θ(-φ) at ω = 0
* Ray reversal is true: θ(θ(φ, pflag, ω), !pflag, ω) ≈ φ
* Magnetic conservation is true: (ω_new(ω_new(ω, pflag), !pflag) ≈ ω
The first property is mandatory and must hold for correct propagation.
The above tests are done for all possible combinations of arguments.

They keyword `only_mandatory` notes whether the rest of
the properties should be tested or not.

## Pinned Particles
In the case of propagation with magnetic field, a particle may be "pinned" (collision-less):
There are no possible collisions that take place and the particle will revolve in circles
forever. This can happen for specific initial conditions depending on your billiard table
and the angular velocity ω.

In such event, the convention followed by `DynamicalBilliards.jl` is the following:
`evolve!()` returns the expected output, however all returned vectors have only 2
entries. The collision times always have the entries `0.0, Inf`. All other returned
vectors have the initial conditions, repeated once.

`evolve!()` can be given an additional `true` optional argument in the case
of magnetic propagation. If you are using ray-splitting it is the 5th argument, otherwise it is the 4th argument. This optional argument throws a `warn()` message whenever a pinned particle is
evolved.

!!! warning "Using `construct()`"
    When using the syntax `construct(evolve!(p, bt, t)...)` be sure that there
    aren't any pinned particles given to evolve. If there are any,
    construct will result in an error.


## Velocity measure

Both `Particle` and `MagneticParticle` are assumed to **always** have a velocity vector of measure 1 during evolution. This simplifies the formulas used internally to a significant amount.

However, during ray-splitting, the a `MagneticParticle` may be in areas with different angular velocities (result of the [ω_new](/tutorials/ray-splitting/#ray-splitting-functions) function). Physically, in such a situation, the velocity measure of the particle could also change. This change depends on the forces acting on the particle (e.g. magnetic field) as well as the relation of the momentum with the velocity (functional type of kinetic energy).

In any case, such a change is not accounted for internally by `DynamicalBilliards`. However it is very easy to implement this by "re-normalizing" the angular velocities you use. Since the "code" velocity has measure one, the rotation radius is given by

```math
r = \frac{1}{\omega_\text{code}} = \frac{v_\text{real}}{\omega_\text{real}}
```

then one simply has to adjust the values of `ω` given in the code with
```math
\omega_\text{code} = \frac{\omega_\text{real}}{v_\text{real}}
```

After getting the timeseries from `construct()`:
```julia
# These are the "code"-data. |v| = 1 always
ct, poss, vels, omegas = evolve!(p, bt, ttotal, ray_splt)
xt, yt, vxt, vyt, t = construct(ct, poss, vels, omegas)
```
you only need to make some final adjustment on the `vxt, vyt`. The position and time data
are completely unaffected.

```julia
omegas_code = omegas
# real angular velocities:
omegas_real = supplied_by_user
# or with some user provided function:
f = o -> (o == 0.5 ? 2o : o*√2)
omegas_real = f.(omegas_code)
# real velocity measure:
vels_real = abs.(omegas_real ./ omegas_code)

contt = cumsum(ct)
omega_t = zeros(t)
vxtreal = copy(vxt)
vytreal = copy(vyt)
j = 1
for i in eachindex(t)
  vxtreal[i] *= vels_real[j]
  vytreal[i] *= vels_real[j]
  omega_t[i] = omegas_real[j]

  if t[i] >=  contt[j]
    j += 1
  end
end
```
Now you can be sure that the particle at time `t[i]` had real velocity `[vxtreal[i], vytreal[i]]` and was propagating with real angular velocity `omega_t[i]`.
