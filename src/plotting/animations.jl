export animate_evolution

# This internal function creates a closure for easy plotting for each particle.
function setup_animation(p::AbstractParticle, bd:: Billiard, t::AbstractFloat,
                         ax::PyPlot.PyCall.PyObject, raysplitters = nothing;
                         dt = 0.01, taillength::Int, pcolor = "k",
                         tailcolor = "C0")

    # get data to plot
    x,y,vx,vy = timeseries(p, bd, t, raysplitters, dt = dt)

    # initial plot
    point = ax[:scatter](x[1], y[1], color = pcolor, s = 30.0, zorder = 2)
    arrow = ax[:quiver](x[1], y[1], 0.08vx[1], 0.08vy[1], color = pcolor,
                        angles = "xy", scale = 1, width = 0.005, zorder = 2)

    tail, = ax[:plot](x[1:2], y[1:2], c = tailcolor, zorder = 1)

    # frame counter
    count = 2
    
    function plot_frame()
        # replot point
        point[:remove]()
        point = ax[:scatter](x[count], y[count], color = pcolor, s = 30.0,
                             zorder = 2)

        # replot arrow
        arrow[:remove]()
        arrow = ax[:quiver](x[count], y[count], 0.08vx[count], 0.08vy[count],
                            color = pcolor, angles = "xy", scale = 1,
                            width = 0.005, zorder = 2)

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
total of `t` times. Optionally enable ray-splitting.

### Keyword Arguments
  * `dt = 0.01` : Time resolution used for production of time series(see 
    [`timeseries`](@ref). It is not recommended to increase this value to 
    preserve the accuracy of the plots.
  * `frameskip = 5` : The amount of `dt`-steps in between succesive frames
  * `tailtime = 1.0` : The length of the "tail" trailing the particle in time 
    units
"""
function animate_evolution(ps::AbstractVector{<:AbstractParticle{T}}, bd::Billiard,
                           t, raysplitters = nothing; dt  = 0.01,
                           frameskip = 5, tailtime = 1.0, ) where {T}

    nps = length(ps)
    taillength = round(Int, tailtime/dt)

    ax = PyPlot.gca()
    plot(bd, ax = ax)
    
    plotframes = Vector{Function}(undef, nps)
    skipframes = Vector{Function}(undef, nps)
    tslengths = Vector{Int}(undef, nps)
    
    for (i, p) ∈ enumerate(ps)
        pf, sf, tl = setup_animation(p, bd, T(t), ax, raysplitters, dt = dt,
                                     taillength = taillength)
        plotframes[i] = pf
        skipframes[i] = sf
        tslengths[i] = tl
    end

    maxframe = minimum(tslengths)

    for i ∈ 2:maxframe
        if (i-1)%frameskip == 0
            [pf() for pf in plotframes]
            sleep(0.01)
        else
            [sf() for sf in skipframes]
        end
    end
end

    
