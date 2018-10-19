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
    tail, = ax[:plot](x[1:2], y[1:2], color = tailcolor; tail_kwargs...)
    
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
total of `t` times. Optionally enable ray-splitting.

### Keyword Arguments
  * `dt = 0.01` : Time resolution used for production of time series(see 
    [`timeseries`](@ref). It is not recommended to increase this value to 
    preserve the accuracy of the plots.
  * `frameskip = 5` : The amount of `dt`-steps in between succesive frames
  * `tailtime = 1.0` : The length of the "tail" trailing the particle in time 
    units
  * `colors` : An array of valid Matplotlib colors for the "tails". If `colors`
    is shorter than `ps`, colors are reused. Defaults to the standard      
    Matplotlib color sequence
  * `particle_kwargs::NamedTuple` : Additional keyword arguments passed to the
    `plot` function for particles.
  * `tail_kwargs::NamedTuple`: Additional keyword arguments passed to the `plot` 
    function for "tails".
"""

function animate_evolution(ps::AbstractVector{<:AbstractParticle{T}},
                           bd::Billiard, t, raysplitters = nothing;
                           dt  = 0.01, frameskip = 5, tailtime = 1.0,
                           colors = nothing, particle_kwargs = NamedTuple(),
                           tail_kwargs = NamedTuple()) where {T}

    nps = length(ps)
    taillength = round(Int, tailtime/dt)

    ax = PyPlot.gca()
    plot(bd, ax = ax)
    
    plotframes = Vector{Function}(undef, nps)
    skipframes = Vector{Function}(undef, nps)
    tslengths = Vector{Int}(undef, nps)

    # choose PyPlot palette by default
    if colors == nothing
        colors = ["C$i" for i ∈ 0:9]
    end
    
    for (i, p) ∈ enumerate(ps)
        
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

    for i ∈ 2:maxframe
        if (i-1)%frameskip == 0
            [pf() for pf in plotframes]
            sleep(0.01)
        else
            [sf() for sf in skipframes]
        end
    end
end

    
