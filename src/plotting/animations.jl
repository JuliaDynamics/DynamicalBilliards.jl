export animate_evolution

# This internal function creates a closure for easy plotting for each particle.
# arguments are mostly similar to animate_evolution
function setup_animation(p::AbstractParticle, bd:: Billiard, t::AbstractFloat,
                         ax::PyPlot.PyCall.PyObject, raysplitters = nothing;
                         dt = 0.01, taillength::Int, tailcolor = "C0",
                         particle_kwargs = NamedTuple(),
                         tail_kwargs = NamedTuple())

    # get data to plot
    x,y,vx,vy = timeseries(p, bd, t, raysplitters, dt = dt)

    # initial plot
    point, arrow = plot_particle(x[1], y[1], vx[1], vy[1]; ax = ax, zorder = 20,
                                 particle_kwargs...)
    tail, = ax[:plot](x[1:2], y[1:2], zorder = 1, color = tailcolor; tail_kwargs...)

    # frame counter
    count = 2

    function plot_frame()
        # replot point
        point[:remove]()
        # replot arrow
        arrow[:remove]()

        point, arrow = plot_particle(x[count], y[count], vx[count], vy[count];
                                     ax = ax, zorder = 20, particle_kwargs...)
        # set tail data
        @views tail[:set_xdata](x[clamp(count-taillength, 1, count):count])
        @views tail[:set_ydata](y[clamp(count-taillength, 1, count):count])

        # increment frame counter
        count += 1
    end

    function skip_frame()
        count += 1
    end

    return plot_frame, skip_frame, length(x)

end

"""
    animate_evolution(ps, bd, colnumber [, raysplitters]; kwargs...)
Animate the evolution of a vector of particles `ps` in billiard `bd` for a
total time `t` (always considered float). Optionally enable ray-splitting.

### Evolution kwargs
  * `dt = 0.01` : Time resolution used for production of time series (see
    [`timeseries!`](@ref). It is not recommended to significantly increase this value,
    to preserve the smoothness of the orbits.
  * `frameskip = 5` : The amount of `dt`-steps performed each frame.
    Increasing either `frameskip` and `dt` makes the animation progress faster.
  * `tailtime = 1.0` : The length of the "tail" trailing the particle in time
    units.
  * `resetting = reset_billiard!` : function called after evolving each individual
    particle in the billiard (so that ray-splitting doesn't brake).
### Colors & plotting kwargs
  * `colors` : An array of valid Matplotlib colors for the "tails". If `colors`
    is shorter than `ps`, colors are reused. Defaults to the standard
    Matplotlib color sequence.
  * `particle_kwargs::NamedTuple` : Additional keyword arguments passed to the
    `plot` function for particles.
  * `tail_kwargs::NamedTuple`: Additional keyword arguments passed to the `plot`
    function for "tails" (line plot).
### Exporting and axis kwargs
  * `figsize = (7.2, 7.2))` : Size for new figure (if one is created).
    Must be divisible by 2 if you want to save the animation.
  * `ax = (figure(figsize = figsize); plot(bd); gca())` : axis to plot on.
  * `savename = nothing` : If given the animation is exported to
    mp4 file (requires ffmpeg). The name can include path.
  * `disable_axis = false` : Remove the axis splines.
  * `deletefigs = true` : To create the animation a lot of figures are saved in
    the save directory and are deleted after the animation is done. You can choose
    to keep them.
  * `dpi = 100` : dpi of saved figures.
  * `framerate = 20` : Animation framerate.
"""
function animate_evolution(ps::AbstractVector{<:AbstractParticle{T}},
                           bd::Billiard, t, raysplitters = nothing;
                           dt  = 0.01, frameskip = 5, tailtime = 1.0,
                           colors = ["C$i" for i ∈ 0:9],
                           particle_kwargs = NamedTuple(),
                           tail_kwargs = NamedTuple(),
                           figsize = (7.2, 7.2),
                           ax = (PyPlot.figure(figsize = figsize); ax = PyPlot.gca(); plot(bd; ax = ax); ax),
                           savename = nothing, dpi = 100, deletefigs = true,
                           disable_axis = false, framerate = 20,
                           resetting = reset_billiard!
                           ) where {T}

    disable_axis && ax[:axis]("off")
    nps = length(ps)
    taillength = round(Int, tailtime/dt)
    savename != nothing && (colnumbers = Int[])

    plotframes = Vector{Function}(undef, nps)
    skipframes = Vector{Function}(undef, nps)
    tslengths = Vector{Int}(undef, nps)

    for (i, p) ∈ enumerate(ps)

        resetting(bd)

        pf, sf, tl = setup_animation(
            p, bd, T(t), ax, raysplitters, dt = dt, taillength = taillength,
            tailcolor = colors[(i-1)%length(colors) + 1],
            tail_kwargs=tail_kwargs, particle_kwargs=particle_kwargs
        )

        plotframes[i] = pf
        skipframes[i] = sf
        tslengths[i] = tl
    end

    maxframe = minimum(tslengths)
    j = 1

    for i ∈ 2:maxframe
        if (i-1)%frameskip == 0
            [pf() for pf in plotframes]
            sleep(0.01)
            if savename != nothing
                s = savename*"_$(j).png"
                PyPlot.savefig(s, dpi = dpi)
                push!(colnumbers, j)
                j += 1
            end
        else
            [sf() for sf in skipframes]
        end

    end

    if savename != nothing
        @assert mod(figsize[1]*dpi, 2) == 0
        @assert mod(figsize[2]*dpi, 2) == 0
        anim = `ffmpeg -y -framerate $(framerate) -start_number 1 -i $(savename)_%d.png
        -c:v libx264 -pix_fmt yuv420p -preset veryslow -profile:v high -level 5.2 $(savename).mp4`
        run(anim)

        if deletefigs
            for i in colnumbers
                rm(savename*"_$(i).png")
            end
        end
    end
end

animate_evolution(p::AbstractParticle, args...; kwargs...) =
animate_evolution([p], args...; kwargs...)
