export bdplot, bdplot!
export bdplot_animstep!
export bdplot_interactive
export bdplot_video
export bdplot_boundarymap

"""
    bdplot(x; kwargs...) → fig, ax
    bdplot!(ax::Axis, x; kwargs...)
Plot an object `x` from `DynamicalBilliards` into a given axis (or a new figure).
`x` can be an obstacle, a particle, a vector of particles, or a billiard.

    bdplot!(ax,::Axis, o::Obstacle; kwargs...)
Keywords are propagated to `lines!` or `poly!`.
Functions `obfill, obcolor, obls, oblw` (not exported)
decide global defaults for linecolor, fillcolor, linestyle, linewidth, when plotting obstacles.

    bdplot!(ax,::Axis, bd::Billiard; clean = true, kwargs...)
If `clean = true`, all axis elements are removed and an equal aspect ratio is establised.
Other keywords are propagated to the obstacle plots.

    bdplot!(ax,::Axis, bd::Billiard, xmin, xmax, ymin, ymax; kwargs...)
This call signature plots periodic billiards: it plots `bd` along its periodic vectors
so that it fills the total amount of space specified by `xmin, xmax, ymin, ymax`.


    bdplot!(ax,::Axis, ps::Vector{<:AbstractParticle}; kwargs...)
Plot particles as a scatter plot (positions) and a quiver plot (velocities).
Keywords `particle_size = 5, velocity_size = 0.05` set the size of plotted particles.
Keyword `colors = JULIADYNAMICS_CMAP` decides the color of the particles, and can be
either a colormap or a vector of colors with equal length to `ps`.
The rest of the keywords are propagated to the scatter plot of the particles.
"""
function bdplot end

"""
    bdplot!(ax, args...; kwargs...)

See [`bdplot`](@ref).
"""
function bdplot! end

"""
    bdplot_boundarymap(bmap, intervals; figkwargs = NamedTuple(), kwargs...)
Plot the output of [`DynamicalBilliards.boundarymap`](@ref) into an axis that
correctly displays information about obstacle arclengths.
Return `figure, axis`.

Also works for the parallelized version of boundary map.

## Keyword Arguments
* `figkwargs = NamedTuple()` keywords propagated to `Figure`.
* `color` : The color to use for the plotted points. Can be either a
  single color or a vector of colors of length `length(bmap)`, in
  order to give each initial condition a different color (for parallelized version).
* All other keywords propagated to `scatter!`.
"""
function bdplot_boundarymap end

"""
    bdplot_interactive(bd::Billiard, ps::Vector{<:AbstractParticle}; kwargs...)
Create a new figure with `bd` plotted, and in it initialize various data for animating the
evolution of `ps` in `bd`. Return `fig, phs, chs`, where `fig` is a figure instance
and `phs, chs` can be used for making custom animations, see below.

## Keywords (interactivity-related)
* `playback_controls = true`: If true, add controls that allow live-interaction with
  the figure, such as pause/play, reset, and creating new particles by clicking
  and dragging on the billiard plot.

## Keywords (visualization-related)
* `dt = 0.001`: The animation always occurs in steps of time `dt`. A slider can decide
  how many steps `dt` to evolve before updating the plots.
* `plot_bmap = false`: If true, add a second plot with the boundary map.
* `colors = JULIADYNAMICS_CMAP` : If a symbol (colormap name) each particle gets
  a color from the map. If Vector of length `N`, each particle gets a color form the vector.
  If Vector with length < `N`, linear interpolation across contained colors is done.
* `tail_length = 1000`: The last `tail_length` positions of each particle are visualized.
* `tail_width = 1`: Width of the dtail plot.
* `fade = true`: Tail color fades away.
* `plot_particles = true`: Besides their tails, also plot the particles as a scatter
  and quiver plot.
* `particle_size = 5`: Marker size for particle scatter plot.
* `velocity_size = 0.05`: Multiplication of particle velocity before plotted as quiver.
* `bmap_size = 4`: Marker size of boundary map scatter plot.
* `figure = NamedTuple()`: Keywords propagated to `Figure` creation.
* `kwargs...`: Remaining keywords are propagated to the billiard plotting.

## Custom Animations
Two helper structures are defined for each particle:
1. `ParticleHelper`: Contains quantities that are updated each `dt` step:
   the particle, time elapsed since last collision, total time ellapsed, tail
   (positions in the last `tail_length` `dt`-sized steps).
2. `CollisionHelper`: Contains quantities that are only updated at collisions:
   index of obstacle to be collided with, time to next collision, total collisions so far,
   boundary map point at upcoming collision.

These two helpers are necessary to transform the simulation into real-time stepping
(1 step = `dt` time), instead of the traditional
DynamicalBilliards.jl setup of discrete time stepping (1 step = 1 collision).

The returned `phs, chs` are two observables, one having vector
of `ParticleHelpers`, the other having vector of `CollisionHelpers`.
Every plotted element is lifted from these observables.

An exported high-level function `bdplot_animstep!(phs, chs, bd, dt; update, intervals)`
progresses the simulation for one `dt` step.
Users should be using `bdplot_animstep!` for custom-made animations,
examples are shown in the documentation online.
The only thing the `update` keyword does is `notify!(phs)`. You can use `false` for it
if you want to step for several `dt` steps before updating plot elements.
Notice that `chs` is always notified when collisions occur irrespectively of `update`.
They keyword `intervals` is `nothing` by default, but if it is `arcintervals(bd)` instead,
then the boundary map field of `chs` is also updated at collisions.
"""
function bdplot_interactive end


"""
    bdplot_video(file::String, bd::Billiard, ps::Vector{<:AbstractParticle}; kwargs...)
Create an animation of `ps` evolving in `bd` and save it into `file`.
This function shares all visualization-related keywords with [`bdplot_interactive`](@ref).
Other keywords are:
* `steps = 10`: How many `dt`-steps are taken between producing a new frame.
* `frames = 1000`: How many frames to produce in total.
* `framerate = 60`.
"""
function bdplot_video end

function bdplot_animstep! end