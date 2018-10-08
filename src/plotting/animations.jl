export animate_evolution

"""
    animate_evolution(p, bd, colnumber [, raysplit]; kwargs...)
Animate the evolution of particle `p` in billiard `bd`, for a total of
`colnumber` collisions. Optionally enable ray-splitting.

### Keyword Arguments
  * `newfig = true` : Creates a new figure at the function call, and plots
    the billiard in that figure.
  * `disable_axis = false` : Turn off the axis of the figure.
  * `framerate = 5` : Rate of drawing new collisions (per second).
    Both during showing animation and during saving it.
  * `col_to_plot = 5` : How many previous collisions are shown during the animation.
  * `particle_kwargs` : A `NamedTuple` of keywords passed to [`plot_particle`](@ref).
  * `orbit_kwargs` : A `NamedTuple` of keywrods passed to `PyPlot.plot`
    which plots the orbit of the particle (`line` object).
  * `savename = nothing` : If `nothing`, nothing is saved.
    If you give a string instead (filename *including path*) an animation `.mp4`
    will be saved with name `savename.mp4`. **This requires
    [`ffmpeg`](https://www.ffmpeg.org) to be accessible from the command line.**
  * `deletefigs = true` : When producing the animation, each frame is saved
    (in the same folder) but then deleted after the animation is created. You
    can choose to keep them by passing `false`.
  * `figsize = (7.2, 7.2)` : `PyPlot` figure size.
  * `dpi = 100` : DPI of saved figures *and* resulting animation.

### Return
The function returns `a, b, c`. Do `a[:remove](); b[:remove](); c[:remove]()` to clear
the particle out of the figure.
"""
function animate_evolution(par::AbstractParticle, bd, colnumber, raysplit = nothing;
    framerate = 5, col_to_plot = 5, savename = nothing,
    particle_kwargs = NamedTuple(), orbit_kwargs = NamedTuple(), newfig = true,
    disable_axis = false, deletefigs = true, dpi = 100,
    figsize = (7.2, 7.2))

    p = copy(par)
    if newfig == true
        fig = PyPlot.figure(figsize = figsize)
        plot(bd; ax = PyPlot.gca())
    end

    disable_axis && PyPlot.axis("off")

    sleeptime = 1/framerate
    i=0
    xdata = Vector{Float64}[]
    ydata = Vector{Float64}[]

    line, point, quiv = nothing, nothing, nothing

    while i < colnumber

        if raysplit == nothing
            xt, yt, vxt, vyt, ts = construct(evolve!(p, bd, 1)...)
        else
            xt, yt, vxt, vyt, ts = construct(evolve!(p, bd, 1, raysplit)...)
        end

        if i < col_to_plot
            push!(xdata, xt)
            push!(ydata, yt)
        else
            popfirst!(xdata); popfirst!(ydata)
            push!(xdata, xt); push!(ydata, yt)
        end

        xpd = Float64[]
        for el in xdata; append!(xpd, el); end

        ypd = Float64[]
        for el in ydata; append!(ypd, el); end

        if i == 0
            line, = PyPlot.plot(xpd, ypd; orbit_kwargs...)
        end
        line[:set_xdata](xpd)
        line[:set_ydata](ypd)

        point, quiv = plot_particle(p; particle_kwargs...)

        if savename != nothing
            s = savename*"_$(i+1).png"
            PyPlot.savefig(s, dpi = dpi)
        end

        sleep(sleeptime)
        if i < colnumber - 1
            point[:remove]()
            quiv[:remove]()
        end
        i+=1
    end
    if savename != nothing
        @assert mod(figsize[1]*dpi, 2) == 0
        @assert mod(figsize[2]*dpi, 2) == 0
        anim = `ffmpeg -y -framerate $(framerate) -start_number 1 -i $(savename)_%d.png
        -c:v libx264 -pix_fmt yuv420p -preset veryslow -profile:v high -level 5.2 $(savename).mp4`
        run(anim)

        if deletefigs
            for i in 1:colnumber
                rm(savename*"_$(i).png")
            end
        end
    end
    return line, point, quiv
end
