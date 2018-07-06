using PyPlot
export animate_evolution

"""
### Keyword Arguments
  * `newfig = true` : Creates a new figure at the function call, and plots
    the billiard in that figure.
  * `disable_axis = false` : Turn off the axis of the figure.
  * `framerate = 5` : Time between frames, both during showing animation and
    during saving it.
  * `col_to_plot = 5` : How many previous collisions are shown during the animation.
  * `particle_kwargs` : Either a Dict{Symbol, Any} or a vector of Tuple{Symbol, Any}.
    Keywords passed into `plot_particle()`.
  * `orbit_kwargs` : Either a Dict{Symbol, Any} or a Vector of Tuple{Symbol, Any}.
    Keywords passed into `PyPlot.plot()` which plots the orbit of the particle
    (`line` object).
  * `savename = nothing` : If `nothing`, nothing is saved.
    If you give a string instead (filename *including path*) an animation `.mp4`
    will be saved with name `savename.mp4`. **This requires the `ffmpeg` to be
    accessible from the command line.**
  * `deletefigs = true` : When producing the animation, each frame is saved
    (in the same folder) but then deleted after the animation is created. You
    can choose to keep them by passing `false`.
  * `figsize = (7.2, 7.2)` : `PyPlot` figure size. Width and height multiplied by
    100 must be divisible by 2.

The function returns `a, b, c`. Do `a[:remove](); b[:remove](); c[:remove]()` to clear
the particle out of the figure.
"""
function animate_evolution(par::AbstractParticle, bd, colnumber, raysplit = nothing;
    framerate = 5, col_to_plot = 5, savename = nothing,
    particle_kwargs = NamedTuple(), orbit_kwargs = NamedTuple(), newfig = true,
    disable_axis = false, deletefigs = true,
    figsize = (7.2, 7.2))

    p = deepcopy(par)
    if newfig == true
        fig = figure(figsize = (7.2, 7.2))
        plot_billiard(bd; ax = gca())
    end

    disable_axis && axis("off")

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
            line, = plot(xpd, ypd; orbit_kwargs...)
        end
        line[:set_xdata](xpd)
        line[:set_ydata](ypd)

        point, quiv = plot_particle(p; particle_kwargs...)

        if savename != nothing
            s = savename*"_$(i+1).png"
            savefig(s, dpi = 100)
        end

        sleep(sleeptime)
        if i < colnumber - 1
            point[:remove]()
            quiv[:remove]()
        end
        i+=1
    end
    if savename != nothing
        anim = `ffmpeg -y -framerate $(framerate) -start_number 1 -i $(savename)_%d.png
        -c:v libx264 -pix_fmt yuv420p -preset veryslow -profile:v baseline -level 3.0 $(savename).mp4`
        run(anim)

        if deletefigs
            for i in 1:colnumber
                rm(savename*"_$(i).png")
            end
        end
    end
    return line, point, quiv
end
