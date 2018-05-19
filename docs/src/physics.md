This page briefly discusses physical aspects of billiard systems.

## Numerical Precision

All core types of `DynamicalBilliards.jl` are parametrically constructed, with
parameter `T<:AbstractFloat`. This means that the fields of all particles and obstacles
contain numbers strictly of type `T`. You will understand why this choice happened
as you continue reading this paragraph.

The main concerns during evolution in a billiard table are:

1. The particle must never leak out of the billiard table. This is simply translated
   to the `distance()` function being **always** positive **after** any collision.
2. The collision time is never infinite, besides the cases of
   [pinned particles](physics/#pinned-particles) in a magnetic billiard.
3. The `relocate!()` process is always finite (and very swift!).

These are solved with the following approach: after the minimum collision time has been calculated, a "test propagation" is done on the position of the particle. If the
`distance(p, obst)` with the colliding obstacle is found to be "wrong", the collision
time is reduced by the `DynamicalBilliards.timeprec(T)` function, with `T` being the parametric type of both the particle and the obstacle. This reduction happens
geometrically, i.e. the time adjusting increment is self-multiplied by 10 at each
recursive relocating step.

* Magnetic propagation and straight propagation both use `timeprec(T)`.
  For both cases the maximum number of geometric iterations is 1 (for the default
  value of `timeprec(T)`.
* This `timeprec` cannot be used with `PeriodicWall` and RaySplitting obstacles with
  `MagneticParticle`s because when relocating forward in time you get progressively
  shallower and shallower angles. This means that you need huge changes in time for even tiny changes in `distance`.
  * For this special case the function `timeprec_forward(T)` is used instead. This
    function results to on average 3-5 geometric relocation steps.

The current definition of these functions is:
```julia
timeprec(::Type{T}) where {T} = eps(T)^(4/5)
timeprec_forward(::Type{T}) where {T} = eps(T)^(3/4)
timeprec_forward(::Type{BigFloat}) = BigFloat(1e-12)
```

Adjusting the global precision of `DynamicalBilliards` is very easy and can be done in
two ways:

1. Choose the floating precision you would like, by initializing your billiard table
   with parametric type `T`, e.g. `bt = billiard_sinai(Float16(0.3))`. This choice
   will propagate to the entire `bt`, all particles resulting from `randominside()`,
   **as well as the entire evolution process**.
2. Re-define the functions `timeprec(T)` and `timeprec_forward(T)`. Decreasing their
   values will make the evolution process slower, but the resulting numbers given by
   `evolve!()` will be more precise. There are no guarantees if you follow this method
   as it may lead to breaking code if `timeprec_forward` becomes too small.

!!! danger "BigFloats"
    Evolution with `BigFloat` in `DynamicalBilliards` is on average
    3 to 4 orders of magnitude slower than with `Float64`.

---


## Pinned Particles
In the case of propagation with magnetic field, a particle may be "pinned" (collision-less):
There are no possible collisions that take place and the particle will revolve in circles
forever. This can happen for specific initial conditions depending on your billiard table
and the angular velocity ω.

In such event, the convention followed by `DynamicalBilliards` is the following:
`evolve!()` returns the expected output, however all returned vectors have only 2
entries. The collision times always have the entries `0.0, Inf`. All other returned
vectors have the initial conditions, repeated once.

`evolve!()` can be given an additional `warning` keyword argument in the case
of magnetic propagation, e.g. `warning = true`. This optional argument throws a `warn()` message whenever a pinned particle is evolved.

!!! warning "Using `construct()`"
    When using the syntax `construct(evolve!(p, bt, t)...)` be sure that there
    aren't any pinned particles given to evolve. If there are any,
    construct will result in an error.

---

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
