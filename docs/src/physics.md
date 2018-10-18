This page briefly discusses physical aspects of billiard systems.

## Pinned Particles
In the case of propagation with magnetic field, a particle may be "pinned" (collision-less):
There are no possible collisions that take place and the particle will revolve in circles
forever. This can happen for specific initial conditions depending on your billiard table
and the angular velocity ω. The function [`ispinned`](@ref) shows you whether a particle meets the conditions.
```@docs
ispinned
```
---

In such event, the convention followed by `DynamicalBilliards` is the following:
[`evolve!`](@ref) returns the expected output, however all returned vectors have only 2
entries. The collision times always have the entries `0.0, Inf`. All other returned
vectors have the initial conditions, repeated once.

[`evolve!`](@ref) can be given an additional `warning` keyword argument in the case
of magnetic propagation, e.g. `warning = true` that throws a  message whenever a pinned particle is evolved.

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

After getting the timeseries:
```julia
# These are the "code"-data. |v| = 1 always
ct, poss, vels, omegas = evolve(p, bd, ttotal)
xt, yt, vxt, vyt, t = timeseries(p, bd, ttotal)
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
