using PyPlot
export plot_boundarymap

"""
    plot_boundarymap(ξs, φs, intervals; kwargs...)

Plots the Poincaré surface of section in boundary coordinates.
The input arguments are the return values of `boundarymap`.

## Keyword Arguments
* `ax = PyPlot.gca()` : The axis to plot on.
* `color = "C0"` : The color to use for the plotted points. Can be either a
  color for `PyPlot.plot` or a vector of colors of length `length(ξs)`, in
  order to give each initial condition a different color.
* `ms = 1.0` : Marker size of the points.
* `bordercolor = "C3"` : The color of the vertical lines that denote the obstacle
  borders.
* `obstacleindices = true`: If true, the a twin axis is created labelling the
  indivdual obstacles by their index
* Any other keyword argument is passed to `PyPlot.plot` which plots the points of
  the section.
"""
function plot_boundarymap(ξs, φs, intervals; ax = PyPlot.gca(),
    color = "C0", bordercolor = "C3", ms = 1.0, obstacleindices = true, kwargs...)

    # Plot PSOS
    for (i, (xis, phis)) in enumerate(zip(ξs, φs))
        c = typeof(color) <: AbstractVector ? color[i] : color

        ax[:plot](xis, phis; marker="o", color = c,
        linestyle="None", ms = ms, kwargs...)
    end

    # Plot obstacle limits
    xmax = 0.0
    for intv ∈ intervals
        for xval ∈ intv
            ax[:plot]([xval,xval], [-π/2, π/2], linewidth = 1.5, color = bordercolor,
            alpha = 0.5)
            xmax = (xval > xmax) ? xval : xmax
        end
    end
    ax[:set_xlim](0,xmax)
    ax[:set_ylim](-π/2,π/2)
    ax[:set_xlabel](L"arc length $\xi$")
    ax[:set_ylabel](L"angle of incidence $\phi$")

    #label tics in units of π/4
    ax[:yaxis][:set_major_formatter](matplotlib[:ticker][:FuncFormatter]((x,p)->"$(x/pi) \$\\pi\$"))
    ax[:yaxis][:set_major_locator](matplotlib[:ticker][:MultipleLocator](base=π/4))

    #number obstacles by index
    if obstacleindices
        #introduce twin axis
        ax2 = ax[:twiny]()

        #non-labelled major tics at every obstacle border
        ax2[:xaxis][:set_major_formatter](matplotlib[:ticker][:NullFormatter]())
        ax2[:set_xticks]([0,map(x->x[2],intervals)...])

        #zero-length minor tics in between borders, labelled with the obstacle index
        ax2[:xaxis][:set_minor_locator](matplotlib[:ticker][:FixedLocator](map(x->mean(x), intervals)))
        ax2[:xaxis][:set_minor_formatter](matplotlib[:ticker][:FixedFormatter](map(x->dec(x),1:length(intervals))))
        ax2[:tick_params](axis="x", which="minor", length=0)#

        ax2[:set_xlabel]("obstacle index")

        return (ax, ax2)
    end

    return ax
end
