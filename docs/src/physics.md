Introductory stuff.

# Implementation

# Ray-Splitting Functions
The condition for transmission is simply: `T(φ, where, ω) > rand()`. If it returns `true`, transmission (i.e. ray-splitting) will happen. As it has already been discussed in the
[Ray-Splitting tutorial](/tutorials/ray-splitting), the condition of total internal reflection
must be taken care by the user.

The [three key functions](/tutorials/ray-splitting/#ray-splitting-functions) given to the `ray-splitter` dictionary must have some properties in order to have physical meaning, like for example that the scattering probability function is even towards φ. One of these properties is absolutely **mandatory** for this package to work properly. This is the property of total internal reflection, i.e. if the refraction angle is calculated to be greater/equal than π/2, no transmission can happen. **This condition is not assured internally** and thefore you must be sure that your transmission probability function satisfies it. In the above example, the function `T` makes sure to return 0 in that case.

In order to test if the `raysplitter` dictionary you have defined has physical meaning, the function `isphysical()` is provided. Its [documentation string](/basic/library/#DynamicalBilliards.isphysical) has all the details one should know.

# Velocity measure

Both `Particle` and `MagneticParticle` are assumed to **always** have a velocity vector of measure 1 during evolution. This simplifies the formulas used internally to a significant amount.

However, during ray-splitting, the a `MagneticParticle` may be in areas with different angular velocities (result of the [ω_new](/tutorials/ray-splitting/#ray-splitting-functions) function). Physically, in such a situation, the velocity measure of the particle could also change. This change depends on the forces acting on the particle (e.g. magnetic field) as well as the relation of the momentum with the velocity (functional type of kinetic energy).

In any case, such a change is not accounted for internally by `DynamicalBilliards`. However it is very easy to implement this by "re-normalizing" the angular velocities you use. Since the "code" velocity has measure one, the rotation radius is given by

`` r = \frac{1}{\omega_\text{code}} = \frac{v_\text{real}}{\omega_\text{real}} ``

then one simply has to adjust the values of `ω` given in the code with ``\omega_\text{code} = \frac{\omega_\text{real}}{v_\text{real}} ``.

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
