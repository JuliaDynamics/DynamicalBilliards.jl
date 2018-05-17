# Ray-Splitting
Ray-splitting is a semi-classical approach to the billiard system, giving a wave attribute to the ray traced by the particle.
Upon collision a particle may propagate through an obstacle (transmission & refraction) or be reflected. Following the mindset of this package, implementing a ray-splitting billiard requires only three very simple steps.

## Ray-Splitting Obstacles
The first step is that an [`Obstacle`](@ref) that supports ray-splitting is required to be present in your billiard table. The only new feature these obstacles have is an additional Boolean field called `pflag` (propagation flag). This field notes on which side of the obstacle the particle is currently propagating.

The normal vector as well as the distance from boundary change sign depending on the value of `pflag`. The obstacles `Antidot` and `SplitterWall` are the equivalents of disk and wall for ray-splitting.

Let's add an `Antidot` to a billiard table:

```julia
using DynamicalBilliards
bt = Obstacle[billiard_rectangle().obstacles...]
a = Antidot([0.5,0.5], 0.3)
push!(bt, a)
bt = Billiard(bt)
```
```
Billiard{Float64} with 5 obstacles:
  Bottom wall
  Right wall
  Top wall
  Left wall
  Antidot
```

## The `RaySplitter` structure
In the second step, you have to define 2+1 functions: transmission probability,
refraction angle and optionally new angular velocity after transmission. These functions, as well as the obstacle index are bundled into a special structure:
```@docs
RaySplitter
```

Notice that if you want different type of transmission/refraction functions for
different obstacles, then you define multiple `RaySplitter`s.

After you have created your `RaySplitter`s, you must bundle them into a tuple
before passing them into a high-level function.

## Ray-Splitter Dictionary
To pass the information of the aforementioned functions into the main API a dictionary is required:
```julia
raysplitter::Dict{Int, Any}
```
This dictionary is a map of the **obstacle index** within the billiard table to a **container of the ray-splitting functions**. This container could be a `Vector` or a `Tuple` and the later is the suggested version.

For example, if we wanted to allocate ray-splitting functions for the **5th** obstacle in our billiard table, which could be e.g. an `Antidot`, we would write something like:
```julia
sa = (θ, pflag, ω) -> pflag ? 2θ : 0.5θ  # refraction (scatter) angle
T = (θ, pflag, ω) -> begin   # Transmission probability
  if pflag
    abs(θ) < π/4 ? 0.5exp(-(θ)^2/2(π/8)^2) : 0.0
  else
    0.75*exp(-(θ)^2/2(π/4)^2)
  end
end
newo = (ω, bool) -> bool ? -0.5ω : -2ω   # new angular velocity
raysplitter = Dict(5 => (T, sa, newo))  # Index maps to container of Functions
```

!!! note "Order of Arguments"
    The functions **must accept the specific number of arguments shown in the previous section** even if some are not used. Also, the functions must be given **in the specific order: (1. transmission probability, 2. refraction angle, 3. new ω)** in the container passed to the dictionary.

The next step is very simple: the `raysplitter` dictionary is directly passed into [`evolve!`](@ref) as a fourth argument.
Using the billiard table we defined previously, where its 5th element is a ray-splitting `Antidot`, we now do:
```julia
ω = 1.25
p = randominside(bt, ω)
dt = 0.05
xt, yt, vxt, vyt, ts = construct(evolve!(p, bt, 25.0, raysplitter)..., dt)
using PyPlot
plot(xt, yt)
```
which should give a result similar to this:

![Ray-splitting billiard](http://i.imgur.com/UfGQfOm.png)

A final difference to be noted is that in the case of ray-splitting with magnetic fields, the fourth value returned by `evolve!()` is not a number, but a vector of angular velocities `omegas`. The value `omegas[i]` is the angular velocity the particle has while propagating from state `pos[i], vel[i]` to state `pos[i+1], vel[i+1]`. The `construct()` function takes this into account.

### No field `pflag` error

If you ever encounter the error `ERROR: type SomeObstacleType has no field pflag` this means that the index provided by your ray-splitting dictionary points to an object that does not support ray-splitting. Use the functions:
```julia
acceptable_raysplitter(raysplitter, bt)
supports_raysplitting(obst::Obstacle)
```
to find out what you did wrong. Most likely, the index you supplied was incorrect, i.e. the index could be `5` instead of `4`.

## Example Animation
In the [examples page](examples), you can find the code for the following animation, which includes ray-splitting:

![Ray-splitter animation](http://i.imgur.com/89s0fon.gif)

## Physics
The condition for transmission is simply: `T(φ, pflag, ω) > rand()`. If it returns `true`, transmission (i.e. ray-splitting) will happen. Otherwise just specular reflection will take place. A more detailed discussion is on the ray-splitting section of the
[Physics](/physics#ray-splitting-functions) page.

## The Ray-Splitting Algorithm
In this section we describe the algorithm we follow to implement the ray-splitting
process. Let $T$ denote the transmission function, $\theta$ the refraction function and $\omega_{\text{new}}$ the new angular velocity function. The following describes the process after a particle has reached an obstacle that supports ray-splitting.

1.  Find the angle of incidence $\phi' = \pi - \arccos(\vec{v} \cdot \vec{n}) = \arccos(\vec{v} \cdot (-\vec{n}))$ with $\vec{n}$ the normal vector at collision point. Notice that we use here $-\vec{n}$ because the velocity is opposing the normal vector before the collision happens. Using $-\vec{n}$ gives the angle between 0 and $\pi/2$ instead of $\pi/2$ to $\pi$.

2.  Find the *correct* sign of the incidence angle, $\phi = \pm \phi'$. Specifically, use the cross product: if the third entry of $\vec{v} \times \vec{n}$ is negative, then have minus sign. The "correct" sign debates on whether the velocity vector is to the right or to the left of $(-\vec{n})$. This is important for finding the correct transmission angle and/or probability.

3.  Check if $T(\phi, \verb|pflag|, \omega) > \text{random}()$. If not, do standard specular reflection.

4. If ray-splitting happens, then relocate the particle so that it is on the *other* side of the colliding obstacle. This contrasts the main evolution algorithm of this billiard package.

5. Re-compute the *correct* angle of incidence, as the position of the particle generally changes with relocating.

6.  Find refraction angle $\theta(\phi, \verb|pflag|, \omega)$. Notice that this is a relative angle with respect to the normal vector. Also notice that $\theta$ may have opposite sign from $\phi$. It depends on the user if they want to add anomalous refraction.

7.  Set `obstacle.pflag = !obstacle.pflag` for *all* obstacles affected by the current `RaySplitter`. This reverses $\vec{n}$ to $-\vec{n}$ as well! So from now on $\vec{n}$ is the opposite than what it was at the beginning of the algorithm!

8.  Find the refraction angle in absolute space. First find $a = \text{atan2}(n_y, n_x)$ and then set $\Theta = a + \theta$.

9. Perform refraction, i.e. set the particle velocity to the direction of $\Theta$.

10. Scale the magnetic field, i.e. set `p.omega` = $\omega_{\text{new}}(\omega, \verb|!pflag|)$. It is important to note that we use `!pflag` because we have already changed the `pflag` field.
